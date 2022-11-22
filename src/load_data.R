x1.start= 3
# shift to start at x1:
# we fix x0 according to first observed angle, which is at index 2 (t0=2)
# first data is produced at x1, at t1=t0+1=3

get_fname = function(key="all"){
  l = list()
  l["lab_data"] =  "mousedata.mat"
  l["par_lab_dir"] = "data_res/"
  l["sim_dir"] = "sim_res/res_2022-11-18-1748/"
  
  if(key=="all")return(names(l))
  if(!key %in% names(l))return()
  
  return(as.character(l[key]))
}

# load result from simulation
load_sim <- function(nst,nst_true){
  library(R.matlab)
  library(stringr)
  
  sim_dir = get_fname("sim_dir")
  in_sim_dir = list.files(sim_dir)
  key = paste("param",nst,"true",nst_true,sep = "_")
  fname = in_sim_dir[str_detect(in_sim_dir,key)]
  if(is.null(fname)) return()
  
  fname = paste(sim_dir,fname,sep="")
  return(readMat(fname))
}

# the dir we store figures from and tables from simulations
sim_figs <-function(){
  sim_dir = get_fname("sim_dir")
  return(paste(sim_dir,"sim_figs/",sep=""))
}

# load store parameter estimate
load_par_est<-function(nst, neurons=59){
  # when more than one estimate stored: load last stored
  library(R.matlab)
  library(stringr)
  
  par_dir = get_fname("par_lab_dir")
  in_par_dir = list.files(par_dir)
  key = paste(nst,"states",neurons,sep = "_")
  fname = sort((in_par_dir[str_detect(in_par_dir,key)]), T)[1]
  
  if(is.na(fname)) return()
  fname = paste(par_dir,fname,sep="")
   
  return(readMat(fname))
}

# load data of neuron activity
load_cell_data<-function(from.x1=FALSE){
  cell_data = load_lab_data()$celldata
  rsums = rowSums(cell_data)
  thrs = 100
  
  # remove boring data
  cell_data = cell_data[rsums>thrs,]
  
  if(from.x1){
    cell_data = cell_data[,x1.start:ncol(cell_data)]
  }
  
  return(cell_data)
  
}

# load store angles
load_ang_data <-function(from.x1=F){
  angs = as.vector(load_lab_data()$resampledAwakeHeadAngleData)
  if(from.x1)angs=angs[x1.start:length(angs)]
  
  return(angs)
}

# support function: cell_data and ang_data is in 'lab_data'
load_lab_data<-function(){
  library(R.matlab)
  fname = get_fname("lab_data")
  return(readMat(fname))
}

# support function: get indexes were angular data is recorded 
get_not_nan<-function(observed_ang=NULL,from.x1=F){
  if(is.null(observed_ang))observed_ang = load_ang_data(from.x1)
  
  not_nan = (1:length(observed_ang))[!is.na(observed_ang)]
  
  return(not_nan)
}
