import scipy
import scipy.stats
import scipy.io
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.animation
plt.ion()

# Whether to animate the state sequence
# Certain graphics backends are somehow very inefficient
# (VSCode interactive with "%matplotlib tk" is fast)
animate = True

# Import data
fname = "45_states_59_neurons_2022-11-14-162252.mat"
# fname = "25_states_59_neurons_2022-11-01-220446.mat"
# fname = "40_states_59_neurons_2022-11-14-140729.mat"
# fname = "60_states_59_neurons_2022-11-14-144627.mat"
mat = scipy.io.loadmat(fname)
tramat = mat["tramat"]
state_sequence = mat["state_sequence"]
angdata = mat["angdata"]
if len(angdata) == 1:
    angdata = angdata[0]

# Prepare weghts and graph stuff
tramat_no_diag = np.copy(tramat)
np.fill_diagonal(tramat_no_diag, 0)
tramat_no_diag_normal = tramat_no_diag / np.max(tramat_no_diag)
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

if animate:
    fig = plt.figure()
    ax = fig.add_subplot(projection = 'polar')
    ax.set_yticklabels("")
    ax.set_theta_zero_location('N')
    title = ax.text(0.5, 0.5, "Timestep: 0", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5})

    line1, = ax.plot([0, meanvalues[int(state_sequence[0])-1]],[0,1], color = 'b', linewidth = 1)
    point1, = ax.plot(angdata[0],1, color='g', marker='o', markersize = 4)
    shade1 = ax.fill_between(np.linspace(meanvalues[int(state_sequence[0])-1]-stdvalues[int(state_sequence[0])-1], meanvalues[int(state_sequence[0])-1]+stdvalues[int(state_sequence[0])-1],10),0,1, alpha=0.3, color='b')
    shade1.set_animated(True)

    legend_elements = [matplotlib.lines.Line2D([0],[0], color='b', lw=1, label="Inferred angle"),
                    matplotlib.lines.Line2D([0],[0], color='r', lw=1, label="No inferred angle"),
                    matplotlib.lines.Line2D([0],[0], color='w', marker='o', markerfacecolor='g', markersize=8, label="Measured angle"),
                    matplotlib.lines.Line2D([0],[0], color='w', marker='o', markerfacecolor='r', markersize=8, label="No measured angle")]
    ax.legend(handles=legend_elements, loc="center")

    def update(t):
        if meanvalues[int(state_sequence[t])-1] == 0:
            line1.set_color('r')
            shade1 = ax.fill_between(np.linspace(line1.get_data()[0][1]-stdvalues[np.where(meanvalues == line1.get_data()[0][1])[0][0]], line1.get_data()[0][1]+stdvalues[np.where(meanvalues == line1.get_data()[0][1])[0][0]],10),0,1, alpha=0.3, color='r')
        else:
            line1.set_data([0, meanvalues[int(state_sequence[t])-1]],[0,1])
            line1.set_color('b')
            shade1 = ax.fill_between(np.linspace(meanvalues[int(state_sequence[t])-1]-stdvalues[int(state_sequence[t])-1], meanvalues[int(state_sequence[t])-1]+stdvalues[int(state_sequence[t])-1],10),0,1, alpha=0.3, color='b')
        shade1.set_animated(True)
        if np.isnan(angdata[t]):
            point1.set_color('r')
        else:
            point1.set_data([angdata[t],1])
            point1.set_color('g')
        title.set_text("Timestep: " + str(t))
        return line1, point1, title, shade1

    frames = np.arange(0,len(state_sequence))

    fig.canvas.draw()
    ani = matplotlib.animation.FuncAnimation(fig, update, frames=frames, blit=True, interval=20)

    plt.show()

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
                        color = (0,0,0,0.9* 0.5 * (
                            tramat_no_diag_normal[inds[from_state]][inds[to_state]]
                            + tramat_no_diag_normal[inds[to_state]][inds[from_state]]
                            + 0.1))
                    )

edge_colors = [state_network[u][v]['color'] for u,v in state_network.edges()]
edge_widths = 25/np.max(weights) * np.array([state_network[u][v]['weight'] for u,v in state_network.edges()]) + 1

# Extract number of states used in the graph
num_used_states = len(list(state_network.nodes))

# Fix the position of the first node (gives a similar graph rotation every time)
#r = 2 * np.pi / (num_used_states * np.sqrt(num_used_states))
#first_node = list(state_network.nodes)[0]
#fixed_pos = {first_node: [r, 2 * r]}

# Find a spring layout and a spectral layout for the graph
pos1 = nx.spring_layout(
    state_network, #pos=fixed_pos, fixed=[first_node],
    iterations=1000, threshold=1e-40)
pos2 = nx.spectral_layout(state_network)
pos3 = nx.spiral_layout(state_network)

# Plot the graph layouts
plt.figure()
nx.draw(state_network, pos1, with_labels=True, edge_color = edge_colors, width = edge_widths)
plt.show()
plt.figure()
nx.draw(state_network, pos2, with_labels=True, edge_color = edge_colors, width = edge_widths)
plt.show()
plt.figure()
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
plt.show(block=True)