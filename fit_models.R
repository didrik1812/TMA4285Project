library(depmixS4)
# ============================
# depmix on 2-dim response
# ===========================
fit_simple_model <- function(y_mat, ns) {
  # number of rows in y_mat is 2 (or we only fit the model to the two first rows)
  # each row is the response from one neuron

  ntimes = ncol(y_mat)
  resp0 = rep(5, ns * 2)
  instart = runif(ns)
  instart = instart / (sum(instart))
  trstart = matrix(data = rep(instart, ns), byrow = TRUE, nrow = ns)

 mod =  depmix(response = list(y_mat[1,] ~ 1, y_mat[2,] ~ 1), # one response-model-formulation for each row in y_mat
                nstates = ns, # number of states
                family = list(poisson(), poisson()),  # each response is poisson (with default link)
               instart = instart, ntimes = ntimes, trstart = trstart, respstart = resp0) # initial values of parameters

  fmod = fit(mod)
  return(fmod)
}

# ============================
# depmix on k-dim response
# ===========================
fit_model <- function(y_mat, ns) {
  # y_mat: matrix of response, column t gives response at time t
  # ns: number of states

  nneur = nrow(y_mat) # number of neurons
  ntimes = ncol(y_mat) # number of observations


  # init values (should probably fit data better):
  # the lambdas are given as a vector. (some) investigation indicate that this correspond to a matrix
  # constructed by row, with each row corresponding to a state
  # (compare summary to response-elements in getpars)
  resp0 = rep(2, ns * nneur)
  instart = runif(ns)
  instart = instart / (sum(instart))
  trstart = matrix(data = rep(instart, ns), byrow = TRUE, nrow = ns)


  # model formulation for vector of response: a list of  'y_mat[i, ] ~ 1'
  # note: formulation is stored as string. we store 'y_mat[1,]~1', 'y_mat[2,]~1' ...
  modform = list()
  for (ry in 1:nneur){
    # paste 'y_mat[ry,] ~ 1'
    modform[[ry]] = formula(paste("y_mat[", ry, ",]~1", sep = ""))
  }

  # response family for vector of poisson response
  fam = rep(list(poisson()), nneur)

  mod <- depmix(modform, nstates = ns, family = fam,
                instart = instart, ntimes = ntimes, trstart = trstart, respstart = resp0)

  fmod = fit(mod)
  return(fmod)
}




# ============================
# MAKE parameters
# ===========================

head_pos_param <- function(num_state, num_neuron) {
  statenames = paste(rep("St", num_state), 1:num_state, sep = "")
  resnames = paste(rep("Res", num_neuron), 1:num_neuron, sep = "")

  p = diag(0.5, num_state)
  for (rown in 0:(num_state - 1)) {
    i = (rown - 1) %% num_state + 1
    j = (rown + 1) %% num_state + 1

    p[rown + 1, i] = p[rown + 1, j] = 0.25
  }

  colnames(p) = rownames(p) = statenames
  rate_state = matrix(data = rep(2 * 1:num_state, num_neuron), nrow = num_neuron, byrow = TRUE)

  noise_col = rnorm(num_neuron, 0, 1)
  rate_neu = matrix(data = rep(noise_col, num_state), ncol = num_state, byrow = FALSE)

  rate_mat = (rate_neu + rate_state)
  rate_mat = abs(rate_mat) # in case we manage to get negative values

  # reorder rate_mate
  colnames(rate_mat) = statenames
  rownames(rate_mat) = resnames
  order_of_lam = names(sort(rate_mat[1, ], decreasing = TRUE))
  rate_mat = rate_mat[, order_of_lam]
  colnames(rate_mat) = statenames

  return(list(lamdas = rate_mat,
              transmat = p))
}

# ====================
#  make sample
# =====================
pois.HMM.generate_sample = function(lambda, Mtrans, StatDist = NULL, n = 500) {
  # n = data length, m = number of states, Mtrans = transition matrix,
  # StatDist = stationary distn
  m = nrow(Mtrans)
  nneur = nrow(lambda)

  # maybe better to force start on state 1, so statDist = c(1,0,....,0) ?
  if (is.null(StatDist)) StatDist = solve(t(diag(m) - Mtrans + 1), rep(1, m))
  mvect = 1:m
  state = numeric(n)

  # get start state
  state[1] = sample(mvect, 1, prob = StatDist)

  # sample hidden chain
  for (i in 2:n)
    state[i] = sample(mvect, 1, prob = Mtrans[state[i - 1], ])

  # sample data for hidden chain:
  y = matrix(data = rpois(n * nneur, lambda = lambda[, state]), ncol = n) #by columns
  return(list(y = y, state = state))
}

# ============================================
# find response matrix, px0 and transistion matrix
# from getpars(fomd)
# ===========================================
unsorted_pars <- function(pars, ns) {
  # find transition matrix, response matrix, and initial prob. from vector of parameters
  print(ns)
  print(length(pars))
  px0 = pars[1:ns] # first ns element are inital prob.
  tramat = matrix(pars[(ns + 1):(ns + ns * ns)], nrow = ns,byrow = TRUE) # next ns*ns elements are transition matrix. (by row)
  resp = matrix(pars[(ns + ns * ns + 1):length(pars)], byrow = TRUE, nrow = ns) # last element are response matrix (by row)

  return(list(resp = resp, tramat = tramat, px0 = px0))
}



# ================================
# RUN DEPMIX ON CONSTRUCTED DATA
# ================================
# pars = head_pos_param(5, 2) # 5 states, 2 neurons
# cell_data = pois.HMM.generate_sample(pars$lamdas, pars$transmat)$y
# fmod = fit_model(cell_data, 5)
# summary(fmod)

# # NOTE: response is on log scale
# # NOTE2: St in fit could correspond to different St in true parameter
# unsorted_pars(getpars(fmod), 5)



# ================================
# RUN DEPMIX ON REAL DATA
# ================================

# Initialize data location and number of states to fit
file_name = "mousedata.mat"
num_states = 35
quick_test = F # If fewer steps and neurons should be used for testing

# Read data
library(R.matlab)
all_data = readMat(file_name)
cell_data = all_data$celldata
angdata = all_data$resampledAwakeHeadAngleData


## Remove neurons that are not very active
thrfrs = rowSums(cell_data)
THR = 100 # perhaps a good threshold
cell_data = cell_data[thrfrs > THR, ]

# Shorten the set to test the code faster
if (quick_test) {
  neurons = 10
  tsteps = 1000
  cell_data = cell_data[1:neurons, 1:tsteps] # test with less data
  angdata = angdata[1:tsteps]
}

# Fit the model and gather results
fmod = fit_model(cell_data, num_states)
summary(fmod)
results = unsorted_pars(getpars(fmod), num_states)
post = posterior(fmod)

# Export results to mat file
writeMat(paste(toString(num_states), "_states_", toString(dim(cell_data)[1]), "_neurons_", format(Sys.time(), "%Y-%m-%d-%H%M%S.mat"), sep = ""), resp = results$resp, tramat = results$tramat, px0 = results$px0, angdata = angdata, post = post, state_sequence = post$state)
