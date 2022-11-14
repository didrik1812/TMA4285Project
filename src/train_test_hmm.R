source("src/my_hmm.R")

train.test.hmm<- setRefClass("train.test.hmm", fields = c("ntest",
                                                          "frac.train",
                                                          "ngap",
                                                          "test.o"),
                           contains =c("my.hmm"))


train.test.hmm$methods(initialize=function(param,y,frac.train=NULL,ngap=NULL){
  .self$config(frac.train,ngap)
  
  end_train = as.integer(.self$frac.train*ncol(y))
  train.y = y[,1:end_train]
  callSuper(param,train.y,run.filter=T)
  
  not_in_test = 1:(end_train +.self$ngap)
  test.y = y[,-(not_in_test)]
  .self$test.o = my.hmm(param=param,y=test.y,run.filter=F)
  
  .self$test.o$px0 = compute.px0()
  .self$test.o$filter() # run filter with correct px0
  
  .self$ntest = test.o$ntimes
  
})
train.test.hmm$methods(config = function(ft,ng){
  default.train = 0.75
  default.gap = 500
  
  if(is.null(ft)) frac.train <<- default.train
  if(is.null(ng)) ngap <<- default.gap
})
train.test.hmm$methods(paste.config = function(){
  confstr = paste("\nFraction of train: ",frac.train,"\nTime steps gap: ",ngap,"\nTime steps test: ",ntest,sep="")
  return(confstr)
})
train.test.hmm$methods(compute.px0=function(){
  lpx0 = solve(t(diag(nst) - Tr + 1), rep(1, nst))
  return(lpx0)
})
train.test.hmm$methods(compute.test.ll=function(){
  test.o$filter()
  return(test.o$ll)
})

train.test.hmm$methods(update.test.param=function(){
 test.o$Tr <<- Tr
 test.o$lambdas <<- lambdas
 test.o$px0 <<- compute.px0()
})

train.test.hmm$methods(optim.param.test=function(max.it){
  '
  method for optimizing train parameter, and find test ll for
  each iteation.
  '
  start_time <- Sys.time()
  print(paste("start ll:", ll))
  
  for(it in 1:max.it){
    start_it <- Sys.time()
    
    cw=tryCatch(
      expr = {
          capture.output(.self$optim.param(1),file=nullfile())
          },
      warning = function(w){ # break loop with warnings
        message("got a warning: ", conditionMessage(w))
        w
        })
    
    if(inherits(cw,"warning")) break
    
    update.test.param()
    compute.test.ll()
    print(paste("done with iteration",it,"in",format(Sys.time()-start_it)))
    print(paste("train ll:",ll))
    print(paste("test ll :",test.o$ll))
    
  }
  print(paste("done with",it,"iterations in",format(Sys.time()-start_time)))
})