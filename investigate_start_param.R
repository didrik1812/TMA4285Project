source("src/plot_lambdas.R")
nst = 60

lc = lamb.w.cases(nst=nst)
plot_lambdas(lc$lam.h.05, lc$lambdas)
plot_lambdas(lc$lam.l.5, lc$lambdas)
plot_lambdas(lc$lam.diff.5, lc$lambdas)


# par(mfrow=c(3,2))
# numbs=1:5
# for(ind in 1:6){
#   plot_lambdas(numbs, lc$lambdas)
#   numbs = numbs+5
# }




