source("src/train_test_hmm.R")
source("src/start_param.R")
source("src/load_data.R")

store_dir = "train_test_res/"
new_res_f = paste(store_dir,"result-",format(Sys.time(), "%Y-%m-%d-%H%M%S"),".txt",sep="")

start_x1 = T
st_v = c(20,25,30,35,40)
res_v = c()
y = load_cell_data(start_x1)

for(nst in st_v){
  print(paste("TRAIN FOR",nst,"STATES"))
  
  param = start_params(nst)
  mod = train.test.hmm(param,y)
  mod$optim.param.test(30) # run to crash or for 30 iterations 
  
  this_dir = paste(store_dir,"nst_",nst,sep="")
  dir.create(this_dir, showWarnings=T)
  
  filename = paste(this_dir,"/testll_",round(mod$test.o$ll),sep="")
  txt_v = c("test ll:",toString(mod$test.o$ll_v),
            "train ll:",toString(mod$ll_v))
  
  
  
  writeLines(txt_v, paste(filename,".txt",sep=""))
  
  # store some more data 
  writeMat(paste(filename,".mat",sep=""),
           resp = mod$lambdas, 
           tramat = mod$Tr, 
           px0 = mod$px0,
           start_x1=start_x1,
           ll_test = mod$test.o$ll_v,
           ll_train = mod$ll_v)
  
  
  new_res_line = paste("nst:",nst,"tes-ll:",mod$test.o$ll)
  cat(new_res_line, file = new_res_f, sep = "\n",append=T)
  
  res_v = c(res_v,mod$test.o$ll)
}
best_line = paste("best:",st_v[which.max(res_v)])
cat(best_line, file = new_res_f, sep = "\n",append=T)
cat(mod$paste.config(),file = new_res_f,sep="\n",append=T)

par(mfrow=c(1,2))
plot(mod$ll_v[-(1:5)],main="loglike train")
plot(mod$test.o$ll_v[-(1:5)],main="loglike test")
