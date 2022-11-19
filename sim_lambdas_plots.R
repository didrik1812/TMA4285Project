source("src/plot_lambdas.R")
lc = lamb.w.cases(nst=40)
store_result = F


stored_40_25 = load_sim(40,25)$resp
stored_25_40 = load_sim(25,40)$resp
# stored_25_25 = load_sim(25,25)$resp

lambdas25 = matrix(start_params(25)$resp,nrow=nrow(lc$lambdas))

if(store_result){
  dir.create(sim_figs())
  png(file=paste(sim_figs(),"lam_40_25.png",sep=""),width=500, height=400)
  par(mar=c(4, 4,3, 1))
}
par(mfrow=c(2,1))
plot_lambdas(lc$lam.diff.5,lambdas25,main="True parameters 25 states",leg.c = -0.2)
plot_lambdas(lc$lam.diff.5,stored_40_25,main="Estimated parameters 40 states (True: 25)",leg.c = -0.2)

if(store_result) dev.off()

if(store_result){
  png(file=paste(sim_figs(),"lam_25_40.png",sep=""),width=500, height=400)
  par(mar=c(4, 4,3, 1))
}
par(mfrow=c(2,1))
plot_lambdas(lc$lam.diff.5,lc$lambdas,main="True parameters 40 states",leg.c = -0.2)
plot_lambdas(lc$lam.diff.5,stored_25_40,main="Estimated parameters 25 states (True: 40)",leg.c = -0.2)

if(store_result) dev.off()
