source("src/train_test_hmm.R")
source("src/start_param.R")
source("src/load_data.R")

store_dir = "train_test_res/"
timetxt = format(Sys.time(), "%Y-%m-%d-%H%M%S")
new_res_f = paste(store_dir, "result-", timetxt, ".txt", sep = "")

start_x1 = TRUE
st_v = 2:100 #c(20,25,30,25,40)
res_v = c()
all = c()
y = load_cell_data(start_x1)

for (nst in st_v){
  print(paste("TRAIN FOR", nst, "STATES"))

  param = start_params(nst)
  mod = train.test.hmm(param, y)
  mod$optim.param.test(30) # run to crash or for 30 iterations

  this_dir = paste(store_dir, "nst_", nst, sep = "")
  dir.create(this_dir, showWarnings = TRUE)

  filename = paste(this_dir, "/testll_", round(mod$test.o$ll), sep = "")
  txt_v = c("test ll:", toString(mod$test.o$ll_v),
            "train ll:", toString(mod$ll_v))
  all = append(all, c(mod$ll_v, mod$test.o$ll_v))

  writeLines(txt_v, paste(filename, ".txt", sep = ""))

  # # store some more data
  # writeMat(paste(filename,".mat",sep=""),
  #          resp = mod$lambdas,
  #          tramat = mod$Tr,
  #          px0 = mod$px0,
  #          start_x1=start_x1,
  #          angdata = load_ang_data(from.x1 = start_x1),
  #          state_sequence = mod$viterbi(),
  #          ll_test = mod$test.o$ll_v,
  #          ll_train = mod$ll_v)

  new_res_line = paste("nst:", nst, "tes-ll:", mod$test.o$ll)
  cat(new_res_line, file = new_res_f, sep = "\n", append = TRUE)

  res_v = c(res_v, mod$test.o$ll)
}
best_line = paste("10 best:", toString(st_v[order(res_v, decreasing = TRUE)[1:10]])) #paste("best:",st_v[which.max(res_v)])
cat(best_line, file = new_res_f, sep = "\n", append = TRUE)
cat(mod$paste.config(), file = new_res_f, sep = "\n", append = TRUE)

writeLines(txt_v, paste(filename, ".txt", sep = ""))

# Store the state and ll data
writeMat(paste(store_dir, "result-", timetxt, ".mat", sep = ""),
         info = mod$paste.config(),
         st_v = st_v,
         res_v = res_v,
         train_test_data = all)

par(mfrow = c(1, 2))
plot(mod$ll_v[-(1:5)], main = "loglike train")
plot(mod$test.o$ll_v[-(1:5)], main = "loglike test")