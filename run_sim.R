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
source("src/print.R")
source("src/compute_store_sim_stats.R")

# ========================================
#   dir name for new results
# ======================================

time_str = format(Sys.time(), "%Y-%m-%d-%H%M")
this_dir = paste("sim_res/res_",time_str,"/",sep="")
in_rapport = paste(this_dir,"sim_figs/",sep="")

dir.create(this_dir, showWarnings=T)
dir.create(in_rapport, showWarnings=T)

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
cases = c("Est St: 25, True st: 25", "Est St: 40, True st: 25","Est St: 25, True st: 40")

meansmat = matrix(0,nrow=3,ncol=2)
rownames(meansmat) = cases
colnames(meansmat) = c("Mean Self-transition","Mean neighbor transition")

ll_mat = matrix(0,nrow=3,ncol=2)
rownames(ll_mat) = cases
colnames(ll_mat) = c("True ll", "Est ll")

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

rownr = 1
ll_mat[rownr,1] = mod$true.ll; ll_mat[rownr,2] = mod$ll
meansmat[rownr,] = estimate_stats(mod$Tr,truep$tramat,rownr,in_rapport)
store_sim_res(mod,nst)
# ========================================
#   b) nst=25, attempt-to-fit-nst = 40
# ========================================
#
# have tried (stored local) simulation of number of states 40, true param and start param as in above. 
# This has not converted after 47 iteration. 
# more interesting: when true number of states is 25, and we try to fit a model with 40 states
# and also: true number of states is 40 and we try to fit a model with 25
startp = start_params(40)
mod$set.params(startp)
mod$ll_v = mod$true.ll # reset ll_v

mod$optim.param(30) # not converted
plot(mod$ll_v)
plot(mod$ll_v[-2])
plot(mod$ll_v[-c(1:10)])

rownr = 2
ll_mat[rownr,1] = mod$true.ll; ll_mat[rownr,2] = mod$ll
meansmat[rownr,] = estimate_stats(mod$Tr,truep$tramat,rownr,in_rapport)
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

rownr = 3
ll_mat[rownr,1] = mod$true.ll; ll_mat[rownr,2] = mod$ll
meansmat[rownr,] = estimate_stats(mod$Tr,truep$tramat,rownr,in_rapport)
store_sim_res(mod,nst)

# ========================================
#   store sum result
# ========================================
ll_mat_s = ll_mat - ll_mat[1,1]
write.csv(ll_mat_s,file=paste(this_dir,"ll_shifted_",ll_mat[1,1],".csv",sep=""))
write.csv(meansmat,file=paste(this_dir,"mean_results.csv",sep=""))

cat(mat.print(ll_mat_s,d1=0),file=paste(in_rapport,"ll_tab.tex",sep=""))
cat(mat.print(meansmat,d1=3),file=paste(in_rapport,"mean_tab.tex",sep=""))

