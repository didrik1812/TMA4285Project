source("src/load_data.R")
#library(circular)


# ==================================
# transition matrix with structure 
# according to hypothesis
# =================================
make_trans_matrix=function(nst,to.self=0.5, to.next=0.25){
  p = diag(to.self,nst) # set probability to stay
  
  # set probability to neighbor states 
  for(rown in 0:(nst-1)){
    i = (rown-1)%%nst +1 
    j = (rown+1)%%nst +1
    
    p[rown+1,i] = p[rown+1,j] = to.next
  } 
  
  remain_prob = 1 -(to.self+2*to.next)
  distant_prob = remain_prob/(nst-3)
  
  p[p==0] = distant_prob 
  
  return(p)
}

# ==============================
# start parameters naive
# ===============================
# init values (should probably fit data better):
# the lambdas are given as a vector. (some) investigation indicate that this correspond to a matrix 
# constructed by row, whit each row corresponding to a state 
# (compare summary to response-elements in getpars) 
naive_start_params<-function(nst,nneur){
  resp0 = runif(nst*nneur,1,5)
  instart = runif(nst)
  instart = instart/(sum(instart))
  
  trstart= matrix(runif(nst*nst),byrow=T,nrow=nst)
  trstart = trstart/rowSums(trstart)
  
  return(list(resp=resp0,
              tramat=trstart,
              px0=instart))
  
}


# ==============================
# start parameters for lab-data
# ===============================
start_params<-function(nst){
  observed_ang = load_ang_data()
  not_nan = get_not_nan(observed_ang)
  
  observed_ang = observed_ang[not_nan]%%(2*pi)
  observed_map = make_observed_map(observed_ang, nst)
  cell_data = load_cell_data()[,not_nan]
  
  smallest_resp = 0.001 # idea: small, but not to hard to "escape" in optimization 
  
  # set lambdas from observed data for angle interval
  resp_param = matrix(data=NaN,nrow =nrow(cell_data),ncol = nst)
  for(st in 1:nst){
    st_cols = cell_data[,observed_map$map[[st]]]
    if(ncol(st_cols) > 5){
      # set lambda to mean for this state (or smallest value)
      lambda_col = rowMeans(st_cols)
      lambda_col[lambda_col<smallest_resp] = smallest_resp
      
      resp_param[,st] = lambda_col 
    }
  }
  
  # interpolate missing values
  if(sum(is.na(resp_param[1,]))>0){
    # interpolate by row
    for(neur in 1:nrow(resp_param)){
      resp_param[neur,] = circular.interpol(resp_param[neur,])
    } 
  }
  
  # set px0
  px0 = rep(0,nst)
  px0[observed_map$x0] = 1
  
  return(list(resp=as.vector(resp_param),
              tramat=make_start_trans(nst),
              px0 =px0))
}


make_start_trans <-function(nst){
  return(matrix(1/nst,nrow=nst,ncol=nst))
}

make_near_hyp_mat <-function(nst){
  to.self = 0.24
  to.next = 0.12
  p = make_trans_matrix(nst,to.self,to.next)
  
  #not_neigh_prob = 1- 
  #p[p==0]=0.05 #  set a x s.t. 0 < x << 0.25 
  #
  #p = p/rowSums(p) # normalize
  
  return(p)
}

make_observed_map <-function(observed_ang,nst){
  inkr = 2*pi/(nst+1)
  start_ang = 1:nst*inkr
  
  all_ind = 1:length(observed_ang)
  observed_map = list()
  
  # find indexes where the angle fit the last state
  last_ind_sm = all_ind[observed_ang<start_ang[1]] # 0< v < 1*inrk
  last_ind_bg = all_ind[observed_ang>=start_ang[nst]] # v > nst*inkr

  observed_map[[nst]]=c(last_ind_sm,last_ind_bg) 
  
  # find x0:
  v0 = observed_ang[1]
  # test for v0 smaller then first start of interval: then x0 is in the last state
  if(v0<start_ang[1]){
    x0 = nst
    } else{
    x0 = which.min(abs(v0-start_ang)) # find closest start of interval
    if(start_ang[x0]>v0) x0=x0-1  # shift when start of interval is bigger than v0  
  }

  # for each state: find indexes where the angle fit the state
  for(st in 1:(nst-1)){
    is_bigger = all_ind[(observed_ang>=start_ang[st])]
    is_smaller = all_ind[(observed_ang<start_ang[(st+1)])]
    
    observed_map[[st]]=intersect(is_smaller,is_bigger) 
  }
  return(list(map=observed_map,x0=x0))
}


circular.interpol <- function(v){
# linear interpolation where value in end-of-vector 
# is assumed to be similar to start-of-vector
  
  start_l = length(v)
  
  # repeat vector 3 times, use part in middle after interpolation
  v = rep(v,3) 
  x=1:length(v)
  xout = x[is.na(v)]
  
  fit_v = approx(x,v,xout = xout) 
  v[xout] = fit_v$y
  
  return(v[(start_l+1):(2*start_l)])
  
  }

# =================================
# parameters for simulation (old)
# ================================

head_pos_param <- function(num_state,num_neuron){
  statenames = paste(rep("St",num_state),1:num_state,sep="")
  resnames = paste(rep("Res",num_neuron),1:num_neuron,sep="")
  
  p = make_trans_matrix(num_state)
  colnames(p)=rownames(p)=statenames
  
  rate_state = matrix(data=rep(2*1:num_state,num_neuron),nrow=num_neuron,byrow = T)
  #sin_trans = sin(pi*(1:num_state)/(num_state+1))
  
  
  noise_col = rnorm(num_neuron,0,1)
  rate_neu = matrix(data=rep(noise_col,num_state),ncol=num_state,byrow = F)
  
  rate_mat = (rate_neu+rate_state)
  rate_mat = abs(rate_mat) # in case we manage to get negative values
  
  # re-order rate_mat
  # colnames(rate_mat) = statenames
  # order_of_lam = names(sort(rate_mat[1,],decreasing=T))
  # rate_mat = rate_mat[,order_of_lam]
  
  rate_mat = matrix(rate_mat,nrow=num_neuron) # in case rate_mate has become a vector
  colnames(rate_mat) = statenames 
  
  if(num_neuron>1)rownames(rate_mat)=resnames
  
  return(list(resp=rate_mat,
              tramat=p))
}