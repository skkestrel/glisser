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

initial_id_to_index = {}
final_id_to_index = {}

with open(sys.argv[1]) as p:
    npl = int(p.readline().strip())
    for i in range(npl * 4): p.readline()
    n = int(p.readline().strip())

    initial = np.zeros((n, 7))
    for i in range(n):
        xyz = list([float(i) for i in p.readline().split()])
        vxyz = list([float(i) for i in p.readline().split()])
        Id = int(p.readline().split()[0])

        initial[i, 0] = xyz[0]
        initial[i, 1] = xyz[1]
        initial[i, 2] = xyz[2]
        initial[i, 3] = vxyz[0]
        initial[i, 4] = vxyz[1]
        initial[i, 5] = vxyz[2]
        initial[i, 6] = Id
        initial_id_to_index[Id] = i

final_to_initial = {}
initial_to_final = {}

with open(sys.argv[2]) as p:
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

        final_id_to_index[pid] = i
        final_to_initial[i] = initial_id_to_index[pid]
        initial_to_final[initial_id_to_index[pid]] = i

flib = np.zeros(initial.shape[0])
flib2 = np.zeros(initial.shape[0])

with open(sys.argv[3]) as p:
    p.readline()
    for line in p:
        flib[final_id_to_index[int(line.split(',')[0])]] = int(line.split(',')[1])
        flib2[final_id_to_index[int(line.split(',')[0])]] = int(line.split(',')[2])

rcenter = 55.46
# 62.63

print("rv2el")

for i in range(initial.shape[0]):
    initial[i, :6] = util2.rv2el(smass, initial[i, :6])
for i in range(final.shape[0]):
    final[i, :6] = util2.rv2el(smass, final[i, :6])

ft = np.zeros(initial.shape[0], dtype=np.bool)


fig, ax = plt.subplots(1, 2)

for i in range(initial.shape[0]):
    ft[i] = int(final[initial_to_final[i]][8]) != 0
ax[0].scatter(initial[ft, 0], initial[ft, 1], s=0.5, c="r", label="non-survivors after 4Gyr")

for i in range(initial.shape[0]):
    ft[i] = int(final[initial_to_final[i]][8]) == 0
ax[0].scatter(initial[ft, 0], initial[ft, 1], s=0.5, c="k", label="survivors after 4Gyr")
ax[0].legend()

for i in range(final.shape[0]):
    ft[i] = int(final[i][8]) == 0 and not flib[i] and not flib2[i]
ax[1].scatter(final[ft, 0], final[ft, 1], s=0.5, c="lime", label="non-librators")

for i in range(final.shape[0]):
    ft[i] = int(final[i][8]) == 0 and (flib[i] or flib2[i])
ax[1].scatter(final[ft, 0], final[ft, 1], s=0.5, c="k", label="librators")


X = np.linspace(rcenter - 1, rcenter + 1, 300)

ax[1].plot(X, 1-30.11/X, c="k")
ax[0].plot(X, 1-30.11/X, c="k")

ax[1].plot(X, 1-20.7/X, c="k", ls="--")
ax[0].plot(X, 1-20.7/X, c="k", ls="--")


ax[1].legend()
ax[0].legend()

ax[0].set_xlim([rcenter - 1, rcenter + 1])
ax[1].set_xlim([rcenter - 1, rcenter + 1])

ax[0].set_ylim([0, 0.7])
ax[1].set_ylim([0, 0.7])

ax[0].set_title("Initial 1Myr librators")
ax[1].set_title("Final 4Gyr survivors")

ax[0].set_xlabel("a (AU)")
ax[1].set_xlabel("a (AU)")
ax[0].set_ylabel("e")
ax[1].set_ylabel("e")

ax[0].tick_params(axis="y",direction="in")
ax[0].tick_params(axis="x",direction="in")
ax[0].set_axisbelow(True)
ax[1].tick_params(axis="y",direction="in")
ax[1].tick_params(axis="x",direction="in")
ax[1].set_axisbelow(True)
ax[0].grid()
ax[1].grid()









fig, ax = plt.subplots(1, 2)

for i in range(final.shape[0]):
    ft[i] = True
ax[0].scatter(initial[ft, 0], initial[ft, 2] / np.pi * 180, s=1, c="k")

for i in range(final.shape[0]):
    ft[i] = int(final[i][8]) == 0 and (flib[i] or flib2[i])
ax[1].scatter(initial[ft, 1], initial[ft, 2] / np.pi * 180,  s=1, c="k")

ax[0].set_xlim([rcenter - 1, rcenter + 1])
ax[1].set_xlim([0,0.7])
# ax[0].set_ylim([0,50])
# ax[1].set_ylim([0,50])
ax[0].set_ylim([0,20])
ax[1].set_ylim([0,20])
ax[0].axhline(y=1,ls="--")
ax[1].axhline(y=1,ls="--")

ax[0].set_ylabel("i_i (deg)")
ax[0].set_xlabel("a_i (AU)")
ax[1].set_xlabel("e_i")

ax[0].tick_params(axis="y",direction="in")
ax[0].tick_params(axis="x",direction="in")
ax[0].set_axisbelow(True)
ax[1].tick_params(axis="y",direction="in")
ax[1].tick_params(axis="x",direction="in")
ax[1].set_axisbelow(True)
ax[0].grid()
ax[1].grid()


plt.figure()
y, edges = np.histogram(initial[:, 1], bins=100, range=[0, 0.8])
centers = (edges[1:] + edges[:-1]) / 2
y = y / y.sum()
plt.plot(centers, y, '-', label="Initial distribution")

for i in range(final.shape[0]):
    ft[i] = int(final[i][8]) == 0 and (flib[i] or flib2[i])
y, edges = np.histogram(final[ft, 1], bins=100, range=[0, 0.8])
y = y / y.sum()
centers = (edges[1:] + edges[:-1]) / 2
plt.plot(centers, y, '-', label="Final distribution")
plt.ylabel("Fraction of total librating particles")
plt.ylim([0, 0.06])
plt.xlim([0, 0.8])
plt.xlabel("e")
plt.legend()


ax = plt.gca()
ax.axvline(1-30.11/rcenter, c="k", ls="-")
ax.axvline(1-20.7/rcenter, c="k", ls="--")

ax.tick_params(axis="y",direction="in")
ax.tick_params(axis="x",direction="in")
ax.set_axisbelow(True)
ax.grid()

plt.show()
