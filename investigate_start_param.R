source("src/plot_lambdas.R")
nst = 60

lc = lamb.w.cases(nst=nst)

# lambdas < 0.5 for all states
plot_lambdas(lc$lam.h.05, lc$lambdas)

# lambdas > 5 for all states
plot_lambdas(lc$lam.l.5, lc$lambdas)

# lambdas with max - min > 5
plot_lambdas(lc$lam.diff.5, lc$lambdas)


# plot all lambdas
par(mfrow=c(3,2))
numbs=1:5
for(ind in 1:6){
  plot_lambdas(numbs, lc$lambdas)
  numbs = numbs+5
}

par(mfrow=c(3,2))
for(ind in 1:6){
  plot_lambdas(numbs, lc$lambdas)
  numbs = numbs+5
  numbs = numbs[numbs<60]
}




