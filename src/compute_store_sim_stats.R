estimate_stats = function(estMat,TrueMat,casenr,save_dir){
  fnames = c("25_true_25","40_true_25","25_true_40")
  self.prob = TrueMat[1,1]
  neigh.prob = TrueMat[1,2]
  not.neigh.prob = TrueMat[1,3]
  
  
  not_neigh = get_not_neigh_el(estMat)
  neigh =  get_neigh_el(estMat)
  self_tr = diag(estMat)
  
  # store means
  meansr = c(mean(self_tr),mean(neigh))
  
  # make plot
  filename = paste(save_dir,"self_",fnames[casenr],".png",sep="")
  png(file=filename,
      width=400, height=200)
  hist(self_tr,xlab="Probability of self-transitions",main="Distribution of self-transitions")
  abline(v=self.prob, col="blue",lwd=3)
  dev.off()
  
  filename = paste(save_dir,"neighbor_",fnames[casenr],".png",sep="")
  png(file=filename,
      width=400, height=200)
  hist(neigh,xlab="Probability of neighbor-transitions",main="Distribution of neighbor-transitions")
  abline(v=neigh.prob, col="blue",lwd=3)
  dev.off()
  
  filename = paste(save_dir,"remain_",fnames[casenr],".png",sep="")
  png(file=filename,
      width=400, height=200)
  hist(not_neigh,xlab="Probability of not-neighbor-transitions",main="Distribution of not-neighbor-transitions")
  abline(v=not.neigh.prob, col="blue",lwd=3)
  dev.off()
  
  return(meansr)
}

get_not_neigh_el = function(estMat){
  lnst = nrow(estMat)
  sc = 1:lnst
  
  estMat[cbind(sc,prev_neigh_ind(lnst))] = NaN
  estMat[cbind(sc,next_neigh_ind(lnst))] = NaN
  diag(estMat) = NaN
  
  return(estMat[!(is.nan(estMat))])
}

get_neigh_el = function(estMat){
  lnst = nrow(estMat)
  sc = 1:lnst
  
  prev_neigh_el = estMat[cbind(sc,prev_neigh_ind(lnst))]
  next_neigh_el = estMat[cbind(sc,next_neigh_ind(lnst))]
  
  return(c(prev_neigh_el,next_neigh_el))
}

prev_neigh_ind = function(lnst) return(c(2:lnst,1))
next_neigh_ind = function(lnst) return(c(lnst,1:(lnst-1))) 