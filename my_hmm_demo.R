library('plot.matrix')
source("src/my_hmm.R") # hmm-class
source("src/start_param.R") # functions for making start parameters
source("src/load_data.R") # functions for loading data


nst = 25
y = load_cell_data(from.x1 = TRUE) # start from first known state
use.homo = T

param = start_params(nst)
if(use.homo) param$transmat=matrix(data=1/nst,ncol=nst,nrow=nst)

mod = my.hmm(param,y) # runs filter in init, takes 6-7 second.
mod$optim.param(1)

# plot transmat
col <- colorRampPalette(c("blue", "white", "red"))
par(oma=c(1,1,1,1.5))
key=list(cex.axis=0.75)
par(oma=c(1,1,1,1.5))
plot(mod$Tr,col=col,key=key)

# continue optim 5 more 
mod$optim.param(5)

# plot likelihood-history (still not converged, it seems)
plot(mod$ll_v[-1])

# plot Tr mat before after
par(oma=c(1,1,1,1.5))
par(mfrow=c(1,2))
plot(mod$Tr,col=col)
plot(param$tramat,col=col)

# plot lambdas before after
par(oma=c(1,1,1,1.5))
par(mfrow=c(1,2))
plot(mod$lambdas,col=col)
plot(matrix(param$lambdas,ncol=ncol(mod$lambdas)),col=col)

# most likely hidden sequence, and post-prob.
gx = mod$viterbi()
s = mod$smoothing()

par(mfrow=c(1,2))
plot(gx[1:200],pch=20,type="o",cex=0.8,main="first entires of most likely hidden sequence")
i = cbind(gx[1:200],1:200)
plot(s[i],pch=20,cex=0.8, main="prop of most likely hidden sequence")

# load stored model
stored.25 = load_par_est(25)
param.stored = list(resp=exp(stored.25$resp),
               tramat=stored.25$tramat,
               px0=as.vector(stored.25$px0))

mod.stored = my.hmm(param.stored, load_cell_data(from.x1 = F)) # load full cell data

# try one more iteration in optim:
mod.stored$optim.param(1)

# looks like: converges to vanishing state
# which produces NaN. the ll is very bad. but interesting.
print(paste("(ll of mod) - (ll of stored)  : ",mod$ll-mod.stored$ll))

# plot lambdas in stored and new model
par(oma=c(1,1,1,1.5))
par(mfrow=c(1,2))
plot(mod$lambdas, main = "new model",col=col)
plot(mod.stored$lambdas, main = "stored model",col=col)

# plot transmat in stored and new model
par(oma=c(1,1,1,1.5))
par(mfrow=c(1,2))
plot(mod$Tr, main = "new model",col=col)
plot(mod.stored$Tr, main = "stored model",col=col)


