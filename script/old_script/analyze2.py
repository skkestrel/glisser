#!/usr/bin/env python3
import matplotlib
# matplotlib.use("Qt5Agg")

import matplotlib.pyplot as plt
import util
import util2
import matplotlib.colors
import matplotlib
import math
import sys
import numpy as np
import matplotlib.style as style
import numpy as np

# style.use('ggplot')

smass = 4 * math.pi * math.pi / math.pow(365.25, 2)

initial, initial_int = util2.read_state(sys.argv[1])
final, final_int = util2.read_state(sys.argv[2])

final_int = final_int[:50000, :]
final = final[:50000, :]

initial_int = initial_int[:final.shape[0], :]
initial = initial[:final.shape[0], :]




#initial = initial[final_int[:, 0].astype(np.int32), :]
#initial_int = initial_int[final_int[:, 0].astype(np.int32), :]

def plot_orbits():
    a = np.linspace(22, 28, 200)
    plt.plot(a, 1 - 19.2 / a)
    plt.text(26, 0.25, "Perihelion uranus")
    plt.plot(a, 30 / a - 1)
    plt.text(26, 0.15, "Aphelion neptune")

# import pdb; pdb.set_trace()

alive = final_int[:, 0] > -1 # always true

final[final[:, 6] == 0, 6] = final[:, 6].max()

# if final.shape[0] < 20001:
if False:
    '''
    plt.scatter(initial[alive, 0], initial[alive, 1], c=final[alive, 6] / 365.25e6, s=1, norm=matplotlib.colors.LogNorm())
    cbar = plt.colorbar()
    cbar.set_label('survival time (MYr)')
    # plt.title('survival time given initial conditions')
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
    ax.grid()
    plt.xlim([23, 27])
    plt.ylim([0.0001, 0.1])
    plt.xlabel('a_i (AU)')
    plt.ylabel('e_i')
    plot_orbits()
    '''

    plt.figure()
    plt.scatter(initial[alive, 0], initial[alive, 2] / 2 / np.pi * 360, c=final[alive, 6] / 365.25e6, s=1, norm=matplotlib.colors.LogNorm())
    cbar = plt.colorbar()
    cbar.set_label('survival time (MYr)')
    #plt.title('survival time given initial conditions')
    plt.xlim([23, 27])
    plt.ylim([0, 0.1])
    plt.xlabel('a_i (AU)')
    plt.ylabel('i_i (deg)')
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
    ax.grid()

    def cFromFlags(flags):
        nflags = []
        for i in range(len(flags)):
            if flags[i] == 0:
                nflags.append("black")
            elif flags[i] == 2:
                nflags.append("red")
            elif flags[i] == 897:
                nflags.append("green")
            elif flags[i] == 1153:
                nflags.append("blue")
            else:
                print(flags[i], "?")
        return nflags

    plt.figure()
    #plt.title('survival time given initial conditions')
    plt.scatter(initial[alive, 2] / math.pi * 180, final[alive, 6] / 365.25e6, s=1, c=cFromFlags(final_int[alive, 1]))
    plt.ylim([0.00001, 4500])
    plt.xlabel('i_i (deg)')
    plt.yscale('log')
    plt.ylabel('survival time (MYr)')
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
    ax.grid()

    plt.figure()
    #plt.title('survival time given initial conditions')
    plt.scatter(initial[alive, 0], final[alive, 6] / 365.25e6, s=1, c=cFromFlags(final_int[alive, 1]))
    plt.ylim([0.00001, 4500])
    plt.xlabel('a (AU)')
    plt.yscale('log')
    plt.ylabel('survival time (MYr)')
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
    ax.grid()


    plt.figure()
    #plt.title('survival time given initial conditions')
    plt.scatter(initial[alive, 1], final[alive, 6] / 365.25e6, s=1, c=cFromFlags(final_int[alive, 1]))
    plt.ylim([0.00001, 4500])
    plt.xlabel('e_i')
    plt.yscale('log')
    plt.ylabel('survival time (MYr)')
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
    ax.grid()

elif final[:, 6].max() > 0:
    print("f")
    '''
    plt.figure()
    #plt.title('survival time given initial conditions')
    util.logHist(plt.gca(), initial[alive, 1], final[alive, 6] / 365.25e6, 1)
    plt.xlabel('e_i')
    plt.yscale('log')
    plt.ylabel('survival time (MYr)')
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
    ax.grid()
    '''

    fig, ax = plt.subplots(2, 1)
    #plt.title('survival time given initial conditions')
    util.dense_scatter(ax[0], initial[alive, 0], final[alive, 6] / 365.25e6, initial[alive, 1], label="e_i", log_y=True, ymin=1, defaultValue=-1)
    ax[0].set_ylim([1, 1e4])
    ax[0].set_xlim([24.2, 26.2])
    ax[0].plot([24.2, 26.2], [4.5e3, 4.5e3], c="r", linewidth=2)
    ax[0].set_xlabel('a_i (AU)')
    ax[0].set_ylabel('survival time (MYr)')
    ax[0].tick_params(axis="y",direction="in")
    ax[0].tick_params(axis="x",direction="in")
    ax[0].set_axisbelow(True)

    #plt.title('survival time given initial conditions')
    util.logHist(ax[1], initial[alive, 0], final[alive, 6] / 365.25e6, 1)
    ax[1].set_xlabel('a_i (AU)')
    ax[1].set_yscale('log')
    ax[1].set_ylabel('survival time (MYr)')
    ax[1].set_ylim([1, 1e4])
    ax[1].set_xlim([24.2, 26.2])
    ax[1].plot([24.2, 26.2], [4.5e3, 4.5e3], c="r", linewidth=2)
    ax[1] = plt.gca()
    ax[1].tick_params(axis="y",direction="in")
    ax[1].tick_params(axis="x",direction="in")
    ax[1].set_axisbelow(True)
    ax[1].grid()

    plt.figure()
    ax = plt.gca()
    #plt.title('survival time given initial conditions')
    util.logHist(ax, initial[alive, 0], final[alive, 6] / 365.25e6, 1)
    ax.set_xlabel('a_i (AU)')
    ax.set_yscale('log')
    ax.set_ylabel('survival time (MYr)')
    ax.set_ylim([1, 1e4])
    ax.set_xlim([24.2, 26.2])
    ax.plot([24.2, 26.2], [4.5e3, 4.5e3], c="r", linewidth=2)
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
    ax.grid()


    plt.figure()
    c = []
    for i in range(final.shape[0]):
        f = int(final_int[i, 1])
        if f == 0:
            c.append("g")
        elif f == 1153:
            c.append("r")
        elif f == 897:
            c.append("b")
        else:
            c.append("k")
    # import pdb; pdb.set_trace()
    plt.scatter(initial[alive, 0], final[alive, 6] / 365.25e6, s=0.1, c=c)
    plt.scatter([], [], s=0.1, c="r", label="Absorbed by Neptune")
    plt.scatter([], [], s=0.1, c="g", label="Alive")
    plt.scatter([], [], s=0.1, c="b", label="Absorbed by Uranus")
    plt.xlabel('a_initial (AU)')
    plt.xlim([24.2, 26.2])
    plt.ylim([1, 4700])
    plt.legend()
    plt.ylabel('survival time (MYr)')
    plt.yscale("log")
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
    ax.grid()



    plt.figure()
    #plt.title('survival time given initial conditions')
    order = list(range(initial.shape[0]))
    order.sort(key=lambda x: final[x, 6])
    util.dense_scatter(plt.gca(), initial[alive, 0], initial[alive, 1], final[alive, 6] / 365.25e6, 150, 150, logBar=True, label="Survival time (MYr)", order=order, upperLimitColor='red', defaultValue=final[alive, 6].min() / 365.25e6)
    plt.xlabel('a_i (AU)')
    plt.ylabel('e_i')
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
    #ax.grid()


    plt.figure()
    plt.title('survival time given initial conditions')
    util.logHist(plt.gca(), initial[alive, 2] / math.pi * 180, final[alive, 6] / 365.25e6, 1)
    plt.xlabel('i (deg)')
    plt.yscale('log')
    plt.ylabel('survival time (MYr)')
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
    ax.grid()

plt.figure()
alive = final_int[:, 1] == 0
alive = np.logical_and(alive, final[:, 0] < 30)

# alive = np.logical_or(alive, final[:, 6] > 365 * 4.5 * 1e9)

#plt.title('surviving particle final states')

# plt.scatter(final[alive, 0], final[alive, 1], s=5, c="blue", label="Final states")
rect = plt.Rectangle([24.2, 0], 2, 0.05, fill=False, zorder=1, linestyle="--", lw="2.5")
ax = plt.gca()
ax.add_patch(rect)
ax.set_xlim([24.1, 26.3])
ax.set_ylim([0, 0.06])
ax.tick_params(axis="y",direction="in")
ax.tick_params(axis="x",direction="in")
ax.set_axisbelow(True)
ax.grid()

if alive.sum() < 1000:
    for i in range(len(alive)):
        if alive[i]:
            pass
            # plt.plot([initial[i, 0], final[i, 0]], [initial[i, 1], final[i, 1]], lw=1, c="pink", zorder=-100)
plt.scatter(initial[alive, 0], initial[alive, 1], s=5, c="red", label="Initial states")

plt.xlabel('a_i (AU)')
plt.legend()
plt.ylabel('e_i')
plot_orbits()


plt.figure()
alive = final_int[:, 1] == 0
#plt.title('surviving particle final states')
plt.scatter(final[alive, 0], final[alive, 2] / math.pi * 180, s=5, c="red", label="Final states")
ax = plt.gca()
ax.tick_params(axis="y",direction="in")
ax.tick_params(axis="x",direction="in")
ax.set_axisbelow(True)
ax.grid()

if alive.sum() < 1000:
    for i in range(len(alive)):
        if alive[i]:
            plt.plot([initial[i, 0], final[i, 0]], [initial[i, 2] / math.pi * 180, final[i, 2] / math.pi * 180], lw=1, c="pink", zorder=-100)
plt.scatter(initial[alive, 0], initial[alive, 2] / math.pi * 180, s=5, c="blue", label="Initial states")

plt.xlabel('a_i (AU)')
plt.legend()
plt.ylabel('i_i (deg)')


plt.show()
