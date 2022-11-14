source("src/simulation_hmm.R")
source("src/start_param.R")
source("src/load_data.R")


nst = 25
ntimes = 15000
startp = start_params(nst)
startp$resp = startp$resp*0.5

truep = start_params(nst)
truep$tramat = make_near_hyp_mat(nst)

mod = simualtion.hmm(truep,startp,ntimes)
mod$optim.param(5)
plot(mod$ll_v[-2])

library('plot.matrix')
col <- colorRampPalette(c("blue", "white", "red"))
par(mfrow=c(1,2))
par(oma=c(1,1,1,1.5))
key=list(cex.axis=0.75)
par(oma=c(1,1,1,1.5))
plot(mod$Tr,col=col,key=key,main="estimated Tr")
plot(truep$tramat,col=col,key=key,main="True Tr")


par(mfrow=c(1,2))
par(oma=c(1,1,1,1.5))
key=list(cex.axis=0.75)
par(oma=c(1,1,1,1.5))
plot(mod$lambdas,col=col,key=key,main="estimated lam")
plot(matrix(truep$resp,ncol=nst),col=col,key=key,main="True lam")
