# =========================
# simulation class
# =========================
source("src/my_hmm.R")

# might create a 'hard.debug' with this kind of data eval
simualtion.hmm<- setRefClass("simualtion.hmm", fields = list(hidden_x="numeric",
                                                         x.counts="matrix",
                                                         true.ll="numeric"), 
                           contains =c("my.hmm"))

simualtion.hmm$methods(initialize=function(paramTrue, paramStart,ntimes){
  dummy.y = paramTrue$resp
  callSuper(paramTrue,dummy.y,run.filter=F) # call super to set parameters
  .self$set.xy(ntimes)
  .self$ntimes = length(hidden_x)
  
  .self$filter() # compute true ll
  .self$true.ll <- ll
  callSuper(paramStart,y,run.filter=T) # call super with start param
})
simualtion.hmm$methods(set.xy=function(n){
  mvect = 1:nst
  hx = numeric(n+1)
  
  hx[1] = sample(mvect, 1, prob = px0)
  for (i in 2:n){
    hx[i] = sample(mvect, 1, prob = Tr[hx[i - 1], ])
  }
  hx = hx[-1] # x0 not in hidden_x
  y <<- matrix(data = rpois(n * nneur, lambda = lambdas[, hx]), ncol = n)
  hidden_x <<- hx
})