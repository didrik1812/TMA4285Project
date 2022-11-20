" 
  script for runing and storing result for simulated data. We look at 
  a) true nst = 25, attepmt-to-fit-nst = 25
  b) true nst = 25, attempt-to-fit-nst = 40
  c) true nst = 40, attempt-to-fit-nst = 25
  
  true parameter for Transistion matrix has structure in line with hypothesis,
  but with probability not-that-close-to-zero for non-neighboors
  
  true parameter for lambdas is the same as used for inital paramters for fitting model
  to lab-data.
  
  for b) and c) we do use the true parameter for lambdas as start
  for a) we start with lambdas*0.5

"

source("src/simulation_hmm.R")
source("src/start_param.R")
# ========================================
#   dir name for new results
# ======================================

time_str = format(Sys.time(), "%Y-%m-%d-%H%M")
this_dir = paste("sim_res/res_",time_str,"/",sep="")


dir.create(this_dir, showWarnings=T)
# ========================================
#   function for this script
# ======================================
store_sim_res = function(mod,true_nst){
  filename = paste(this_dir,"param_",mod$nst,"_true_",nst,"_",round(mod$ll),sep="")
  writeMat(paste(filename,".mat",sep=""),
           true_nst = true_nst,
           resp = mod$lambdas,
           tramat = mod$Tr,
           px0 = mod$px0,
           ll_v = mod$ll_v,
           true_ll= mod$true.ll,
           ll = mod$ll)
  
}

true_param_sim = function(nst){
  truep = start_params(nst)
  truep$tramat = make_near_hyp_mat(nst)
  return(truep)
}

# ========================================
#  global-ish variables
# ========================================
ntimes = 15000
# ========================================
#   a) nst=25, attempt-to-fit-nst = 25
# =======================================
nst = 25
truep = true_param_sim(nst)
startp = start_params(nst)
startp$resp = startp$resp*0.5

mod = simualtion.hmm(truep,startp,ntimes)
mod$optim.param(25) # could stop by 20, a bit clearer convergence by 25
plot(mod$ll_v[-2])
plot(mod$ll_v[-c(1:10)])

library('plot.matrix')
col <- colorRampPalette(c("blue", "white", "red"))
par(mfrow=c(1,2))
key=list(cex.axis=0.75)
par(oma=c(1,1,1,1.5))
plot(mod$Tr,col=col,key=key,main="estimated Tr")
plot(truep$tramat,col=col,key=key,main="True Tr")


par(mfrow=c(1,2))
par(oma=c(1,1,1,1.5))
key=list(cex.axis=0.75)
plot(mod$lambdas,col=col,key=key,main="estimated lam")
plot(matrix(truep$resp,ncol=nst),col=col,key=key,main="True lam")

store_sim_res(mod,nst)
# ========================================
#   b) nst=25, attempt-to-fit-nst = 40
# ========================================
startp = start_params(40)
mod$set.params(startp)
mod$ll_v = mod$true.ll # reset ll_v

mod$optim.param(30) # not converted
plot(mod$ll_v)
plot(mod$ll_v[-2])
plot(mod$ll_v[-c(1:10)])

store_sim_res(mod,nst)

# ========================================
#   c) nst=40, attempt-to-fit-nst = 25
# ========================================
nst = 40
truep = true_param_sim(nst)
startp = start_params(25)
mod = simualtion.hmm(truep,startp,ntimes)
mod$optim.param(30)

# we see: 25 states explain quite a bit less of the variation
plot(mod$ll_v[-2])
plot(mod$ll_v[-c(2:10)])

store_sim_res(mod,nst)