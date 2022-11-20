estimate_stats = function(estMat,TrueMat,casenr,save_dir){
  save_res = save_dir != nullfile()
  
  fnames = c("25_true_25","40_true_25","25_true_40")
  
  # xlabs:
  xlab.self = "Probability of self-transitions"
  xlab.n ="Probability of neighbor-transitions"
  xlab.nn="Probability of not-neighbor-transitions"
  
  # Mains:
  main.self = "Distribution of self-transitions"
  main.n = "Distribution of neighbor-transitions"
  main.nn = "Distribution of not-neighbor-transitions"
  
  # True probs
  self.prob = TrueMat[1,1]
  neigh.prob = TrueMat[1,2]
  not.neigh.prob = TrueMat[1,3]
  
  
  not_neigh = get_not_neigh_el(estMat)
  neigh =  get_neigh_el(estMat)
  self_tr = diag(estMat)
  
  # means
  meansr = c(mean(self_tr),mean(neigh))
  
  # make plot
  if(save_res){
    filename = paste(save_dir,"self_",fnames[casenr],".png",sep="")
    png(file=filename, width=400, height=200)    
  }
  hist(self_tr,xlab=xlab.self,main=main.self)
  abline(v=self.prob, col="blue",lwd=3)
  if(save_res) dev.off()
  
  if(save_res){
    filename = paste(save_dir,"neighbor_",fnames[casenr],".png",sep="")
    png(file=filename, width=400, height=200) 
  }
  hist(neigh,xlab=xlab.n,main=main.n)
  abline(v=neigh.prob, col="blue",lwd=3)
  if(save_res) dev.off()
  
  if(save_res){
    filename = paste(save_dir,"remain_",fnames[casenr],".png",sep="")
    png(file=filename, width=400, height=200)
  }
  hist(not_neigh,xlab=xlab.nn,main=main.nn)
  abline(v=not.neigh.prob, col="blue",lwd=3)
  if(save_res) dev.off()
  
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