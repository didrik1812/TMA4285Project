x1.start= 3
# shift to start at x1:
# we fix x0 according to first observed angle, which is at index 2 (t0=2)
# first data is produced at x1, at t1=t0+1=3

get_fname = function(key="all"){
  l = list()
  l["lab_data"] =  "mousedata.mat"
  l["par_lab_25"] = "25_states_59_neurons_2022-11-01-220446.mat"
  l["par_lab_35"] = "35_states_59_neurons_2022-11-06-111051.mat"
  l["sim_dir"] = "sim_res/res_2022-11-18-1748/"
  
  if(key=="all")return(names(l))
  if(!key %in% names(l))return()
  
  return(as.character(l[key]))
}

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

sim_figs <-function(){
  sim_dir = get_fname("sim_dir")
  return(paste(sim_dir,"sim_figs/",sep=""))
}

load_par_est<-function(nst){
  library(R.matlab)
  key = paste("par_lab",nst,sep = "_")
  fname = get_fname(key)
  
  if(is.null(fname)) return()
  
  return(readMat(fname))
}

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

load_ang_data <-function(from.x1=F){
  angs = as.vector(load_lab_data()$resampledAwakeHeadAngleData)
  if(from.x1)angs=angs[x1.start:length(angs)]
  
  return(angs)
}

load_lab_data<-function(){
  library(R.matlab)
  fname = get_fname("lab_data")
  return(readMat(fname))
}

get_not_nan<-function(observed_ang=NULL,from.x1=F){
  if(is.null(observed_ang))observed_ang = load_ang_data(from.x1)
  
  not_nan = (1:length(observed_ang))[!is.na(observed_ang)]
  
  return(not_nan)
}
