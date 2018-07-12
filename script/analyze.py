#!/usr/bin/env python3
import matplotlib
matplotlib.use("Qt5Agg")

import matplotlib.pyplot as plt
import util
import matplotlib.colors
import matplotlib
import math
import sys
import numpy as np
import matplotlib.style as style
import numpy as np

style.use('ggplot')

smass = 4 * math.pi * math.pi / math.pow(365.25, 2)

with open(sys.argv[2]) as p:
	npl = int(p.readline().strip())
	for i in range(npl * 4): p.readline()
	n = int(p.readline().strip())

	initial = np.zeros((n, 6))
	for i in range(n):
		xyz = list([float(i) for i in p.readline().split()])
		vxyz = list([float(i) for i in p.readline().split()])
		p.readline()

		initial[i, 0] = xyz[0]
		initial[i, 1] = xyz[1]
		initial[i, 2] = xyz[2]
		initial[i, 3] = vxyz[0]
		initial[i, 4] = vxyz[1]
		initial[i, 5] = vxyz[2]

with open(sys.argv[1]) as p:
	npl = int(p.readline().strip())
	for i in range(npl * 4): p.readline()
	n = int(p.readline().strip())

	final = np.zeros((n, 9))
	for i in range(n):
		xyz = list([float(i) for i in p.readline().split()])
		vxyz = list([float(i) for i in p.readline().split()])
		flags = p.readline().strip().split()

		dtime = float(flags[2])
		pid = int(flags[0])

		final[i, 0] = xyz[0]
		final[i, 1] = xyz[1]
		final[i, 2] = xyz[2]
		final[i, 3] = vxyz[0]
		final[i, 4] = vxyz[1]
		final[i, 5] = vxyz[2]
		final[i, 6] = dtime
		final[i, 7] = pid
		final[i, 8] = int(flags[1])
	initial = initial[final[:, 7].astype(np.int32), :]

print("rv2el")

for i in range(initial.shape[0]):
	initial[i, :6] = util.rv2el(smass, initial[i, :6])
for i in range(final.shape[0]):
	final[i, :6] = util.rv2el(smass, final[i, :6])

def plot_orbits():
    a = np.linspace(22, 28, 200)
    plt.plot(a, 1 - 19.2 / a)
    plt.text(26, 0.25, "Perihelion uranus")
    plt.plot(a, 30 / a - 1)
    plt.text(26, 0.15, "Aphelion neptune")


# import pdb; pdb.set_trace()

alive = final[:, 8] > -1
final[final[:, 6] == 0, 6] = final[:, 6].max()

if final.shape[0] < 20001:
    plt.scatter(initial[alive, 0], initial[alive, 1], c=final[alive, 6] / 365.25e6, s=1, norm=matplotlib.colors.LogNorm())
    cbar = plt.colorbar()
    cbar.set_label('survival time (MYr)')
    plt.title('survival time given initial conditions')
    plt.xlim([23, 27])
    plt.ylim([0.0001, 0.1])
    plt.xlabel('a (AU)')
    plt.ylabel('e')
    plot_orbits()

    plt.figure()
    plt.scatter(initial[alive, 0], initial[alive, 2] / 2 / np.pi * 360, c=final[alive, 6] / 365.25e6, s=1, norm=matplotlib.colors.LogNorm())
    cbar = plt.colorbar()
    cbar.set_label('survival time (MYr)')
    plt.title('survival time given initial conditions')
    plt.xlim([23, 27])
    plt.ylim([0, 0.1])
    plt.xlabel('a (AU)')
    plt.ylabel('i (deg)')

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
    plt.title('survival time given initial conditions')
    plt.scatter(initial[alive, 0], final[alive, 6] / 365.25e6, s=1, c=cFromFlags(final[alive, 8]))
    plt.ylim([0.00001, 4500])
    plt.xlabel('a (AU)')
    plt.yscale('log')
    plt.ylabel('survival time (MYr)')

    plt.figure()
    plt.title('survival time given initial conditions')
    plt.scatter(initial[alive, 2] / math.pi * 180, final[alive, 6] / 365.25e6, s=1, c=cFromFlags(final[alive, 8]))
    plt.ylim([0.00001, 4500])
    plt.xlabel('i (deg)')
    plt.yscale('log')
    plt.ylabel('survival time (MYr)')

    plt.figure()
    plt.title('survival time given initial conditions')
    plt.scatter(initial[alive, 1], final[alive, 6] / 365.25e6, s=1, c=cFromFlags(final[alive, 8]))
    plt.ylim([0.00001, 4500])
    plt.xlabel('e')
    plt.yscale('log')
    plt.ylabel('survival time (MYr)')

elif final[:, 8].max() > 0:
    plt.figure()
    plt.title('survival time given initial conditions')
    util.logHist(initial[alive, 1], final[alive, 6] / 365.25e6, 1)
    plt.xlabel('e')
    plt.yscale('log')
    plt.ylabel('survival time (MYr)')

    plt.figure()
    plt.title('survival time given initial conditions')
    util.logHist(initial[alive, 0], final[alive, 6] / 365.25e6, 1)
    plt.xlabel('a (AU)')
    plt.yscale('log')
    plt.ylabel('survival time (MYr)')

    plt.figure()
    plt.title('survival time given initial conditions')
    util.logHist(initial[alive, 2] / math.pi * 180, final[alive, 6] / 365.25e6, 1)
    plt.xlabel('i (deg)')
    plt.yscale('log')
    plt.ylabel('survival time (MYr)')

plt.figure()
alive = final[:, 8] == 0
plt.title('surviving particle final states')
plt.scatter(final[alive, 0], final[alive, 1], s=5, c="red", label="Initial states")

if alive.sum() < 1000:
    for i in range(len(alive)):
        if alive[i]:
            plt.plot([initial[i, 0], final[i, 0]], [initial[i, 1], final[i, 1]], lw=1, c="pink", zorder=-100)
plt.scatter(initial[alive, 0], initial[alive, 1], s=5, c="blue", label="Final states")

plt.xlabel('a (AU)')
plt.legend()
plt.ylabel('e')
plot_orbits()


plt.figure()
alive = final[:, 8] == 0
plt.title('surviving particle final states')
plt.scatter(final[alive, 0], final[alive, 2] / math.pi * 180, s=5, c="red", label="Initial states")

if alive.sum() < 1000:
    for i in range(len(alive)):
        if alive[i]:
            plt.plot([initial[i, 0], final[i, 0]], [initial[i, 2] / math.pi * 180, final[i, 1] / math.pi * 180], lw=1, c="pink", zorder=-100)
plt.scatter(initial[alive, 0], initial[alive, 2] / math.pi * 180, s=5, c="blue", label="Final states")

plt.xlabel('a (AU)')
plt.legend()
plt.ylabel('i (deg)')


plt.show()
