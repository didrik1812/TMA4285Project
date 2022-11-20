source("src/plot_lambdas.R")
source("src/print.R")
source("src/compute_store_sim_stats.R")

# cases: "25_true_25","40_true_25","25_true_40"

lc = lamb.w.cases(nst=40)
store_result = T
save_dir = nullfile() # for estimate_stats function
this_dir = get_fname("sim_dir")

if(store_result){
  save_dir = sim_figs()
  dir.create(sim_figs())
}

# ====================================================
#   TRUE AND ESTIMATED PARAMETERS
# ====================================================
stored_40_25 = load_sim(40,25)
stored_25_40 = load_sim(25,40)
stored_25_25 = load_sim(25,25)


lambdas25 = matrix(start_params(25)$resp,nrow=nrow(lc$lambdas))
Tr25 = make_near_hyp_mat(25)
Tr40 = make_near_hyp_mat(40)

# ====================================================
#  TABLES
# ====================================================  
base_ll = stored_25_25$true.ll # to shift ll with

meansmat = matrix(0,nrow=3,ncol=2)
ll_mat = matrix(0,nrow=3,ncol=2)

# row names
cases = c("Est St: 25, True St: 25", 
          "Est St: 40, True St: 25",
          "Est St: 25, True St: 40")

rownames(meansmat) = cases
rownames(ll_mat) = cases

# col names mean table
colnames(meansmat) = c("Mean Self-transition",
                       "Mean Neighbor-transition")

# col names ll table
colnames(ll_mat) = c("True ll",
                     "Est ll")

# ====================================================
# FILL IN TABELS AND MAKE TRANSIITION STATS FIGURES
# ====================================================
for(casenr in 1:3){
  stored_o = list(stored_25_25,stored_25_40,stored_40_25)[[casenr]]
  true_Tr = list(Tr25,Tr25,Tr40)[[casenr]]
  
  ll_mat[casenr,] = c(stored_o$true.ll,stored_o$ll)
  meansmat[casenr,] =  estimate_stats(stored_o$tramat,true_Tr,casenr,save_dir)
}
rm(stored_o, true_Tr)
# ========================================
#   # STORE TABLES
# ========================================
if(store_result){
  ll_mat_s = ll_mat - base_ll
  #write.csv(ll_mat_s,file=paste(this_dir,"ll_shifted_",base_ll,".csv",sep=""))
  #write.csv(meansmat,file=paste(this_dir,"mean_results.csv",sep=""))
  
  cat(mat.print(ll_mat_s,d1=0),file=paste(sim_figs(),"ll_tab.tex",sep="")) # shifted ll
  cat(mat.print(ll_mat,d1=0),file=paste(sim_figs(),"ll_tab_native.tex",sep="")) # native ll
  cat(mat.print(meansmat,d1=3),file=paste(sim_figs(),"mean_tab.tex",sep="")) # mean stats
}

# ========================================
#  LAMBDAS FIGURES
# ========================================

# fig titles
main.true.25 = "True parameters 25 states"
main.est.40 = "Estimated parameters 40 states (True: 25)"

main.true.40 = "True parameters 40 states"
main.est.25 = "Estimated parameters 25 states (True: 40)"


if(store_result){
  png(file=paste(sim_figs(),"lam_40_25.png",sep=""),width=500, height=400)
  par(mar=c(4, 4,3, 1))
}
par(mfrow=c(2,1))
plot_lambdas(lc$lam.diff.5,lambdas25,main=main.true.25, leg.c = -0.2)
plot_lambdas(lc$lam.diff.5,stored_40_25$resp, main=main.est.40, leg.c = -0.2)

if(store_result) dev.off()

if(store_result){
  png(file=paste(sim_figs(),"lam_25_40.png",sep=""),width=500, height=400)
  par(mar=c(4, 4,3, 1))
}
par(mfrow=c(2,1))
plot_lambdas(lc$lam.diff.5,lc$lambdas, main=main.true.40, leg.c = -0.2)
plot_lambdas(lc$lam.diff.5,stored_25_40$resp, main=main.est.25, leg.c = -0.2)

if(store_result) dev.off()
