source("src/start_param.R")
plot_lambdas<-function(lamnr,inputlamb,leg.c=-0.08,...){
  lambs = inputlamb[lamnr,]
  nneur = nrow(lambs)
  clrs = rainbow(nneur)
  ylim = c(min(lambs),max(lambs))
  
  plot(lambs[1,],type="o",col=clrs[1],pch=16,ylim=ylim,ylab=expression(lambda),xlab="state",...)
  for(neur in 2:nneur){
    lines(lambs[neur,],type="o",col=clrs[neur],pch=16)
  }
  lamnr = as.numeric(lamnr)
  
  l=list()
  for(i in 1:length(lamnr)){
    neur = lamnr[i]
    l[[i]] =bquote(Lambda[.(neur)])
  }
  
  legend("topright", legend=l,pch=16,col=clrs,horiz=T,xpd=TRUE,inset=c(0,leg.c),bty="n",seg.len=1)
  
}

lamb.w.cases = function(lambdas=NULL,nst=40){
  return_full = F
  
  if(is.null(lambdas)) {
    lambdas = matrix(start_params(nst)$resp,ncol=nst)
    return_full = T
  }
  all.lamb = 1:nrow(lambdas)
  l = list()
  
  l[["lam.h.05"]] = all.lamb[array_max(lambdas)<0.5]
  l[["lam.diff.5"]] = all.lamb[array_diff(lambdas)>5]
  l[["lam.l.5"]]=all.lamb[array_min(lambdas)>5]
  
  if(return_full) l[["lambdas"]] = lambdas
  
  return(l)
}

array_min = function(m,bydim=1){
  apply(m,bydim,min)
} 

array_max = function(m,bydim=1){
  apply(m,bydim,max)
} 
array_std = function(m,bydim=1){
  apply(m,bydim,sd)
}
array_range = function(m,bydim=1){
  apply(m,bydim,range)
}
array_diff = function(m,bydim=1){
  v_diff =  function(v){return(max(v)-min(v))}
  apply(m,bydim,v_diff)
}



# lam.h.10 = all.lamb[array_max(lambdas)<10]
# lam.h.5 = all.lamb[array_max(lambdas)<5]
# lam.h.2 = all.lamb[array_max(lambdas)<2]



#plot_lambdas(lam.h.10)
#plot_lambdas(lam.h.5)
#plot_lambdas(lam.h.2)