# TMA4285 Time Series Project: Uncovering Neural Population States in Rat Data
In this project we have used a hidden markov model (HMM) to analyze neural popluation states in rat data. In this repository lies all code used for producing the results of our report. 

## Main code:
Following below are a description of the main files used in our project.

**Plotting**

Here follows the files that are used to produce the figures displayed in the report. 
* `plot_states.py`: plots the log-likelihood against the number of states using results from `run_train_test.R`. I used for selection of number of states.

* `plot_results.py`: Plots spring-graph layout of the states, transition matrix and possibly an animation of the predicted angle vs the actual angle.

* `plot_multiple_results.py`: Plotting spring-graphs and transition matrices for 4 different cases of number of states, coloring of edges is added to symbolizes the strength/weight of the edge.

**Data Fitting**

* `run_train_test.R`: Trains HMM model on data for a different number of states. Optimization is done with respect to the (log)likelihood. Both the loglikelihood of testset and of the trainset is written to a .mat file.
* `run_sim.R`: Trains and test models on simulated model. Results are plotted and written to a .mat file.
* `generate_sim_results.R`: More testing of model through simulation, results are shown in the appendix.
* `investigate_start_param.R`: Plot some figures displaying the effect of different inital conditions.
* `my_hmm_demo.R`: Exploration of implemented HMM model.
* `my_hmm_test.R`: More testing and exploration of implemented HMM model


**Data**

* `mousedata.mat`: Real rat data used in model fitting.


## Other folders

### train_test_res

This folder contains results from model training and testing.

### data_res

This folder contains data from model fitting and is used for producing figures displayed in the report. 

### figs
Here lies figures produced from the code and used in the report.

### src
This folder contains the implementation of the hidden markov model.
Files:
* `my_hmm.R`: Implementation of our HMM through optimization and methods such as filtering, smoothing and the vitberi alogortihm. 
* `train_test_hmm.R`: Builds upon `my_hmm.R` and fits model to a trainset and compares with a testset.
* `load_data.R`: Extract data from .mat files and prep it for model fitting.



## Predicted angle vs actual angle

Animation of infered angle and actual angle for 45 states.

[Animation](https://user-images.githubusercontent.com/92478930/203111350-fb41de37-fb72-40e3-8450-cf9501848908.mp4)


Confidence interval is shown through shaded area.





