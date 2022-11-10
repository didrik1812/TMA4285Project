import scipy
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Import data
fname = "35_states_59_neurons_2022-11-06-111051.mat"
mat = scipy.io.loadmat(fname)
tramat = mat["tramat"]
state_sequence = mat["state_sequence"]
angdata = mat["angdata"][0]

# Prepare weghts and graph stuff
tramat_no_diag = np.copy(tramat)
np.fill_diagonal(tramat_no_diag, 0)
weights = tramat / np.mean(tramat_no_diag)
num_states = len(tramat)
meanvalues = np.zeros(num_states)
stdvalues = np.zeros(num_states)
counts = np.zeros(num_states)
state_network = nx.Graph()

# Find the average head angle per state
for i in range(num_states):
    count = np.sum(state_sequence == i + 1)
    angs = angdata[state_sequence == i + 1]
    if sum(np.isnan(angs)) > 0:
        angs = angs[~np.isnan(angs)]

    # Require at least 10 angle values to calculate average and std
    if len(angs) > 10:
        meanvalues[i] = scipy.stats.circmean(angs)
        stdvalues[i] = scipy.stats.circstd(angs)
    else:
        print(
            "Fewer than 10 of presorted state "
            + str(i + 1)
            + " have an angle associated. No average angle assigned."
        )
    counts[i] = count

# Sort the states by average angle (and counts if tied (only realistically occurs if the mean angle is not set))
inds = np.argsort(
    np.array(
        list(zip(meanvalues, counts)), dtype=[("meanvalues", "f8"), ("counts", "i4")]
    ),
    order=["meanvalues", "counts"],
)
meanvalues = meanvalues[inds]
stdvalues = stdvalues[inds]
counts = counts[inds]

# Print state number, mean angle, angle std, and number of occurences
for i in range(len(inds)):
    print(i + 1, meanvalues[i], stdvalues[i], counts[i])

# Create a transition matrix sorted by average angle
sortedTM = np.zeros(np.shape(tramat))
sortedTM[:] = np.nan
for i in range(len(inds)):
    if counts[i] < 10:
        continue
    for j in range(len(inds)):
        if counts[j] < 10:
            continue
        sortedTM[i, j] = tramat[inds[i], inds[j]]

# Add weighted edges to the graph object
for from_state in range(len(tramat)):
    if counts[from_state] > 10:
        for to_state in range(len(tramat)):
            if to_state != from_state and counts[to_state] > 10:
                if (from_state + 1, to_state + 1) not in state_network.edges:
                    state_network.add_edge(
                        from_state + 1,
                        to_state + 1,
                        weight=0.5
                        * (
                            weights[inds[from_state]][inds[to_state]]
                            + weights[inds[to_state]][inds[from_state]]
                        ),
                    )

# Extract number of states used in the graph
num_used_states = len(list(state_network.nodes))

# Fix the position of the first node (gives a similar graph rotation every time)
r = 2 * np.pi / (num_used_states * np.sqrt(num_used_states))
first_node = list(state_network.nodes)[0]
fixed_pos = {first_node: [r, 2 * r]}

# Find a spring layout and a spectral layout for the graph
pos1 = nx.spring_layout(
    state_network, pos=fixed_pos, fixed=[first_node], iterations=200, threshold=1e-25
)
pos2 = nx.spectral_layout(state_network)
pos3 = nx.spiral_layout(state_network)

# Plot the two graph layouts
nx.draw(state_network, pos1, with_labels=True)
plt.show()
nx.draw(state_network, pos2, with_labels=True)
plt.show()
nx.draw(state_network, pos3, with_labels=True)
plt.show()


# Plot the sorted transition matrix
plt.figure()
im = plt.imshow(
    sortedTM[-num_used_states:, -num_used_states:],
    vmin=0,
    vmax=1,
    interpolation="None",
    cmap="Greys",
    extent=[
        num_states - num_used_states + 0.5,
        num_states + 0.5,
        num_states + 0.5,
        num_states - num_used_states + 0.5,
    ],
)
t = np.linspace(0, 1, num=6)
plt.colorbar(im, ticks=t)
plt.title("Transition matrix")
plt.xlabel("To state")
plt.ylabel("From state")
plt.show()

# Plot the sorted transition matrix with a logarithmic colorbar
plt.figure()
im = plt.imshow(
    sortedTM[-num_used_states:, -num_used_states:],
    norm=LogNorm(vmin=1e-8, vmax=1),
    interpolation="None",
    cmap="Greys",
    extent=[
        num_states - num_used_states + 0.5,
        num_states + 0.5,
        num_states + 0.5,
        num_states - num_used_states + 0.5,
    ],
)
t = np.logspace(-8, 0, num=9)
plt.colorbar(im, ticks=t)
plt.title("Transition matrix (logarithmic colorbar)")
plt.xlabel("To state")
plt.ylabel("From state")
plt.show()
