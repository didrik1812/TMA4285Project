import numpy as np
import scipy
import matplotlib.pyplot as plt

# Load data
train_test_results = scipy.io.loadmat("./train_test_res/result-2022-11-16-103935.mat")
info = train_test_results["info"]
states = train_test_results["st_v"]
test_results = train_test_results["res_v"]
all_data = train_test_results["train_test_data"]

# The all_data got concatenated :/
# Here's a hacky way to separate it again
train_ll = []
test_ll = []
tmp_train = [0]
for e in all_data:
    if len(train_ll) == len(test_ll):
        if tmp_train[0] == 0:
            tmp_train[0] = e
            continue
        if abs((e-tmp_train[-1]) / e) < 0.1:
            tmp_train.append(e)
            continue
        else:
            train_ll.append(tmp_train)
            tmp_test = [e]
            continue
    else:
        if e == all_data[-1] and len(train_ll) == len(states):
            tmp_test.append(e)
            test_ll.append(tmp_test)
            continue
        if abs((e-tmp_test[-1]) / e) < 0.1:
            tmp_test.append(e)
            continue
        else:
            test_ll.append(tmp_test)
            tmp_train = [e]
            continue

# Find the best ll result per state number (not just the final ll)
# It seems the optimization on the train data sometimes worsens 
# the test data results
test_results_best = [max(test_ll[i]) for i in range(len(states))]

plt.figure()
plt.title(info[0])
plt.xlabel("Number of states fitted")
plt.ylabel("loglikelihood on test data")
plt.plot(states, test_results, label = "Final ll")
plt.plot(states, test_results_best, label = "Maximum ll")
plt.legend()
plt.tight_layout()
plt.show()