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
# Saving the animation further requires ffmpeg or imagemagick
# and does some update frame tricks that are much slower than
# running the code with save_anim = False. In addition, all frames
# of this slower execution must be rendered before the code continues. 
# WARNING: For some reasons the color updates don't work properly
# on the displayed figure if save_anim = True, but the actually
# saved animation is correct. ¯\_(ツ)_/¯
animate = False
save_anim = False 

# Import data for 4 cases (Very ungeneral)
fname1 = "15_states_59_neurons_2022-11-16-165359.mat"
fname2 = "30_states_59_neurons_2022-11-16-165712.mat"
fname3 = "45_states_59_neurons_2022-11-16-170235.mat"
fname4 = "60_states_59_neurons_2022-11-16-171208.mat"
mats = [scipy.io.loadmat(fname1),
        scipy.io.loadmat(fname2),
        scipy.io.loadmat(fname3),
        scipy.io.loadmat(fname4)]

tramats = []
state_sequences = []
graphs = []
sortedTMs = []
positions = []
widths = []
colors = []
anims = []

for model in range(4):
    tramats.append(mats[model]["tramat"])
    state_sequences.append(mats[model]["state_sequence"])
    angdata = mats[model]["angdata"]
    if len(angdata) == 1:
        angdata = angdata[0]

    # Prepare weghts and graph stuff
    tramat_no_diag = np.copy(tramats[model])
    np.fill_diagonal(tramat_no_diag, 0)
    tramat_no_diag_normal = tramat_no_diag / np.max(tramat_no_diag)
    weights = tramats[model] / np.mean(tramat_no_diag)
    num_states = len(tramats[model])
    meanvalues = np.zeros(num_states)
    stdvalues = np.zeros(num_states)
    counts = np.zeros(num_states)
    graphs.append(nx.Graph())

    # Find the average head angle per state
    for i in range(num_states):
        count = np.sum(state_sequences[model] == i + 1)
        angs = angdata[state_sequences[model] == i + 1]
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

        line1, = ax.plot([0, meanvalues[int(state_sequences[model][0])-1]],[0,1], color = 'b', linewidth = 1)
        point1, = ax.plot(angdata[0],1, color='g', marker='o', markersize = 4)
        shade1 = ax.fill_between(np.linspace(meanvalues[int(state_sequences[model][0])-1]-stdvalues[int(state_sequences[model][0])-1], meanvalues[int(state_sequences[model][0])-1]+stdvalues[int(state_sequences[model][0])-1],10),0,1, alpha=0.3, color='b')
        shade1.set_animated(True)

        legend_elements = [matplotlib.lines.Line2D([0],[0], color='b', lw=1, label="Inferred angle"),
                        matplotlib.lines.Line2D([0],[0], color='r', lw=1, label="No inferred angle"),
                        matplotlib.lines.Line2D([0],[0], color='w', marker='o', markerfacecolor='g', markersize=8, label="Measured angle"),
                        matplotlib.lines.Line2D([0],[0], color='w', marker='o', markerfacecolor='r', markersize=8, label="No measured angle")]
        ax.legend(handles=legend_elements, loc="center")

        def update(t):
            if save_anim:
                ax.collections.pop()
            if meanvalues[int(state_sequences[model][t])-1] == 0:
                line1.set_color('r')
                shade1 = ax.fill_between(np.linspace(line1.get_data()[0][1]-stdvalues[np.where(meanvalues == line1.get_data()[0][1])[0][0]], line1.get_data()[0][1]+stdvalues[np.where(meanvalues == line1.get_data()[0][1])[0][0]],10),0,1, alpha=0.3, color='r')
            else:
                line1.set_data([0, meanvalues[int(state_sequences[model][t])-1]],[0,1])
                line1.set_color('b')
                shade1 = ax.fill_between(np.linspace(meanvalues[int(state_sequences[model][t])-1]-stdvalues[int(state_sequences[model][t])-1], meanvalues[int(state_sequences[model][t])-1]+stdvalues[int(state_sequences[model][t])-1],10),0,1, alpha=0.3, color='b')
            shade1.set_animated(True)
            if np.isnan(angdata[t]):
                point1.set_color('r')
            else:
                point1.set_data([angdata[t],1])
                point1.set_color('g')
            title.set_text("Timestep: " + str(t))
            return line1, point1, title, shade1

        frames = np.arange(0,len(state_sequences[model]))

        fig.canvas.draw()
        if save_anim:
            ani = matplotlib.animation.FuncAnimation(fig, update, frames=frames[:200], blit=True, interval=50)
            ani.save('animation.svg')#, writer='imagemagick')#, fps=30)
        else:
            ani = matplotlib.animation.FuncAnimation(fig, update, frames=frames, blit=True, interval=20)
            plt.show()
        anims.append(ani)
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
    sortedTM = np.zeros(np.shape(tramats[model]))
    sortedTM[:] = np.nan
    for i in range(len(inds)):
        if counts[i] < 10:
            continue
        for j in range(len(inds)):
            if counts[j] < 10:
                continue
            sortedTM[i, j] = tramats[model][inds[i], inds[j]]
    sortedTMs.append(sortedTM)

    # Add weighted edges to the graph object
    for from_state in range(len(tramats[model])):
        if counts[from_state] > 10:
            for to_state in range(len(tramats[model])):
                if to_state != from_state and counts[to_state] > 10:
                    if (from_state + 1, to_state + 1) not in graphs[model].edges:
                        graphs[model].add_edge(
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

    edge_widths = 25/np.max(weights) * np.array([graphs[model][u][v]['weight'] for u,v in graphs[model].edges()]) + 1
    widths.append(edge_widths)
    edge_colors = [graphs[model][u][v]['color'] for u,v in graphs[model].edges()]
    colors.append(edge_colors)

    # Find a spring layout for the graph
    pos = nx.spring_layout(graphs[model], iterations=1000, threshold=1e-40)
    positions.append(pos)

# Plot the graph layouts
fig, ax = plt.subplots(2, 2)
for i in range(len(graphs)):
    nstates = len(tramats[i])
    ix = np.unravel_index(i, ax.shape)
    plt.sca(ax[ix])
    nx.draw(graphs[i], positions[i], with_labels=True, node_size = 150, font_size = 10, edge_color = colors[i], width = widths[i])
    ax[ix].set_title(str(nstates) + " states", fontsize=10, pad=0)
    ax[ix].set_axis_off()
plt.subplots_adjust(wspace=0,
                    hspace=0.1)
plt.show()

# Plot the transition matrices
fig, ax = plt.subplots(2, 2)
for i in range(len(sortedTMs)):
    nstates = len(tramats[i])
    ix = np.unravel_index(i, ax.shape)
    plt.sca(ax[ix])
    # draw all nodes homogeneously, and edge weights as filtered
    num_used_states = len(list(graphs[i].nodes))
    im = plt.imshow(
        sortedTMs[i][-num_used_states:, -num_used_states:],
        vmin=0,
        vmax=1,
        interpolation="None",
        cmap="Greys",
        extent=[
            nstates - num_used_states + 0.5,
            nstates + 0.5,
            nstates + 0.5,
            nstates - num_used_states + 0.5,
        ],
    )
    ax[ix].set_title(str(nstates) + " states", fontsize=10, pad=1)
    ax[ix].tick_params(axis='both', labelsize=10)
    ax[ix].yaxis.get_major_locator().set_params(integer=True)
    ax[ix].xaxis.set_major_locator(plt.MaxNLocator(8))
    ax[ix].yaxis.set_major_locator(plt.MaxNLocator(8))
fig.supxlabel("To state", fontsize = 10)
fig.supylabel("From state", fontsize = 10)
# fig.suptitle("Transition matrices")
t = np.linspace(0, 1, num=6)
fig.colorbar(im, ticks=t, ax = ax.ravel().tolist())
# plt.subplots_adjust(hspace = 1)
# plt.tight_layout()
plt.show()

# Plot the transition matrices with a logarthmic colorbar
fig, ax = plt.subplots(2, 2)
for i in range(len(sortedTMs)):
    nstates = len(tramats[i])
    ix = np.unravel_index(i, ax.shape)
    plt.sca(ax[ix])
    # draw all nodes homogeneously, and edge weights as filtered
    num_used_states = len(list(graphs[i].nodes))
    im = plt.imshow(
        sortedTMs[i][-num_used_states:, -num_used_states:],
        norm=LogNorm(vmin=1e-8, vmax=1),
        interpolation="None",
        cmap="Greys",
        extent=[
            nstates - num_used_states + 0.5,
            nstates + 0.5,
            nstates + 0.5,
            nstates - num_used_states + 0.5,
        ],
    )
    ax[ix].set_title(str(nstates) + " states", fontsize=10, pad=1)
    ax[ix].tick_params(axis='both', labelsize=10)
    ax[ix].yaxis.get_major_locator().set_params(integer=True)
    ax[ix].xaxis.set_major_locator(plt.MaxNLocator(8))
    ax[ix].yaxis.set_major_locator(plt.MaxNLocator(8))
fig.supxlabel("To state", fontsize = 10)
fig.supylabel("From state", fontsize = 10)
# fig.suptitle("Transition matrices")
t = np.logspace(-8, 0, num=9)
fig.colorbar(im, ticks=t, ax = ax.ravel().tolist())
# plt.subplots_adjust(hspace = 1)
# plt.tight_layout()
plt.show()