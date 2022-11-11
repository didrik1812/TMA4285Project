my.hmm <- setRefClass("my.hmm",field= list(f="matrix",
                                           p="matrix",
                                           s="matrix",
                                           global.x = "numeric",
                                           y ="matrix",
                                           fy = "matrix",
                                           ll="numeric",
                                           ll_v = "numeric",
                                           nst = "numeric",
                                           nneur = "numeric",
                                           ntimes = "numeric",
                                           Tr="matrix",
                                           px0="numeric",
                                           lambdas="matrix",
                                           Tr.hat="matrix",
                                           lambdas.hat = "matrix",
                                           Tr.c="numeric",
                                           Tr.cc = "numeric" ,
                                           Tr.debug = "matrix",
                                           Tr.change="numeric"))
my.hmm$methods(initialize=function(param,y){
  .self$Tr <- param$transmat
  .self$px0 = param$px0
  .self$nst = nrow(Tr)
  .self$lambdas <-matrix(param$lambdas,ncol=nst) # lambdas as matrix 
  .self$nneur = nrow(lambdas) 
  .self$y = matrix(y,nrow=nrow(lambdas)) # data as matrix
  .self$ntimes = ncol(y)
  .self$filter() # run filter 
})

#====================
# optim param
# ===================
my.hmm$methods(optim.param=function(niter){
  start_time <- Sys.time()
  print(paste("start ll:",ll))
  for(it in 1:niter){
    prev.ll = ll
    .self$smoothing()
    Tr.change <<- c(Tr.change,sum(abs(Tr - Tr.hat)))
    
    Tr.copy = Tr 
    lambdas.copy = lambdas

    .self$update.param()
    .self$filter()
    
    if(ll < prev.ll){
      print("divergens: use prev estimate")
      Tr <<- Tr.copy
      lambdas <<- lambdas.copy
      print("rerun filter:")
      .self$filter()
      break
    }
    
    print("")
    print(paste("iteration ",it," ll: ",ll,sep=""))
  }
  end_time= Sys.time()
  print(paste("done with", it,"iterations in",format(end_time-start_time)))
})

my.hmm$methods(update.param=function(){
  Tr <<- Tr.hat
  lambdas <<- lambdas.hat
})
# =========================================== 
#         VITERBI
#===========================================
my.hmm$methods(viterbi=function(){
  max.x = matrix(0,ncol=ntimes,nrow=nst)
  c.indx = 1:nst
  
  log.T = log(Tr)
  log.fy = log(fy)
  
  # forward recursion: find arg.max mu
  mu = log(f[,1])
  for(t in 2:ntimes){
    mu.Tr = mu+log.T
    max.t = apply(mu.Tr,2,which.max) # arg. for highest col values
    
    mu = mu.Tr[cbind(max.t,c.indx)] # correct col.values from mu.Tr
    mu = log.fy[,t]+mu
    #mu = mu/sum(mu)
    
    max.x[,(t-1)] = max.t
  }
  # backward: picking out the correct max.t
  gx = numeric(ntimes)
  xt.max = which.max(mu)
  for(t in (ntimes-1):1){
    gx[t+1]=xt.max
    xt.max = max.x[xt.max, t] # correct value for x_t, is in row xt.max
  }
  gx[1] = xt.max # xt.max is now x1.max
  global.x <<- gx
  return(gx)
  
})


# =========================================== 
#         LESS STUPID SMOOTH
#===========================================
my.hmm$methods(smoothing = function(){
  # NOTE phi.next, is  k*f(y_{t+2:n}|x_t) => phi = k*f(y_{t+1:n}|x_t)
  # (due to normalization)
  #
  # it is still possible for phi -> 0 :
  # since we normalize, its how f(y_{t+1:n}|x_t = i) compares to the other states that matters,
  #
  # f(y_{t+1:n}|x_t = i) =  f(x_t = i|y_{t+1:n})f(y_{t+1:n})/p(x_t = i)
  #                 prop.to f(x_t = i|y_{t+1:n})/p(x_t = i)
  #
  # for this to -> 0 we must (almost) be able to conclude that x_t can't be i, 
  # from the evidence of state t+1, t+2, ...,n
  # Thus the transition probabilities from i, can only be high to states which
  # the evidence find really unlikely to be the true state. Thus: 
  # f(x_t = i,x_{t+1}=j|y_{t:n}) -> 0 for all j, so  it should be ok to set NaN to 0 ?
  #
  
  start_time <- Sys.time()
  phi.next = rep(1,nst)
  ls = matrix(ncol=ntimes,nrow = nst)
  ls[,ntimes] = f[,ntimes]
  Tr.est = matrix(nrow=nst,ncol=nst,data=0)

  for(t in (ntimes-1):1){
    T.fy.phi = t(t(Tr)*fy[,t+1]*phi.next) # '*' goes by columns
    phi = rowSums(T.fy.phi)    
    
    mi = phi*f[,t]
    mi = mi/sum(mi)
    
    Tr.est = Tr.est + T.fy.phi*mi/phi  # k*phi.next/k*phi = phi.next/phi
    Tr.est[is.nan(Tr.est)] = 0  # NaN should mean 0
    
    ls[,t] = mi
    phi.next = phi/sum(phi) # normalize phi, so it does not vanish
  }
  ls[is.nan(ls)] = 0 # NaN should mean 0
  
  s <<- ls
  Tr.debug <<- Tr.est
  Tr.sum = rowSums(Tr.est)
  Tr.hat <<- Tr.est/Tr.sum
  
  # CONTROL
  # Row sums in transition matrix should equal 1 
  # Row sums in Tr.est should equal row sums of smoothing: 
  #  sum_k( sum_t( p(x_t=i, x_{t+1} = k |y) ) ) = sum_t( p(x_t=i|y) )
  s.sums = rowSums(ls[,-ntimes])
  Tr.c <<- c(Tr.c,sum(abs(1-rowSums(Tr.hat)))) 
  Tr.cc <<- c(Tr.cc,sum(abs(Tr.sum - s.sums))) 
  
  s.sums = s.sums + ls[,ntimes] # full sum
  
  # estimate lambdas
  l.lambdas = matrix(0,nrow=nneur,ncol=nst)
  for(neur in 1:nneur){
    y.s.prod = t(s)*y[neur,]
    l.lambdas[neur,] = colSums(y.s.prod)/s.sums
  }
  lambdas.hat <<- l.lambdas
  print(paste("done with smoothing in",format(Sys.time()-start_time)))
  return(ls)
})

# =========================================== 
#         FILTER
#===========================================
my.hmm$methods(filter=function(store.fy=T){
    fy <<- .self$compute_eval()
   
    start_time <- Sys.time()
    
    lT = t(Tr)
    lf = lp = matrix(0,nrow=nst,ncol=(ntimes+1))
    lp[,1] = lT%*%px0
    ll <<- 0
     
    
    # p_j = sum_{i}( Tr_{ij}*f_{i} )
    #     = sum_{i}( t(Tr)_{ji})*f_{i})
    #     = (t(Tr)f)_{j}
    
    for(tp in 1:ntimes){
      ft = lp[,tp]*fy[,tp]
      
      k = sum(ft)
      lf[,tp] = ft/k
      ll <<- ll+log(k)
      
      lp[,(tp+1)] = lT%*%lf[,tp] # will be computed one more time than needed.
    }
    
    f <<- lf[,1:ntimes]
    p <<- lp[,1:ntimes]
    ll_v <<- c(ll_v,ll)
    
    if(!store.fy) fy <<- matrix(ncol = 0,nrow=0)  # fy will vanish from memory when function end?
    print(paste("computed f and p in",format(Sys.time()-start_time)))
   
     return(f)
})

# =========================================== 
#         COMPUTE EVAL 
#===========================================    
my.hmm$methods(compute_eval=function(){
  start_time <- Sys.time()
  
  ntimes <<-ncol(y)
  local_m <- matrix(0,nrow=nst,ncol=ntimes)
  
  
  for(st in 1:nst){
    full_state = dpois(y,lambdas[,st])
    local_m[st,] = apply(full_state,2,prod)
  }
  print(paste("done with eval in",format(Sys.time()-start_time)))
  
  return(local_m)
})