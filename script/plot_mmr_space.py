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

lib = np.zeros(initial.shape[0])
lib2 = np.zeros(initial.shape[0])
flib = np.zeros(initial.shape[0])
flib2 = np.zeros(initial.shape[0])

with open(sys.argv[3]) as p:
    p.readline()
    for line in p:
        lib[initial_id_to_index[int(line.split(',')[0])]] = int(line.split(',')[1])
        lib2[initial_id_to_index[int(line.split(',')[0])]] = int(line.split(',')[2])
with open(sys.argv[4]) as p:
    p.readline()
    for line in p:
        flib[final_id_to_index[int(line.split(',')[0])]] = int(line.split(',')[1])
        flib2[final_id_to_index[int(line.split(',')[0])]] = int(line.split(',')[2])

print("rv2el")

for i in range(initial.shape[0]):
    initial[i, :6] = util.rv2el(smass, initial[i, :6])
for i in range(final.shape[0]):
    final[i, :6] = util.rv2el(smass, final[i, :6])

ft = np.zeros(initial.shape[0], dtype=np.bool)

fig, ax = plt.subplots(1, 2)
for i in range(initial.shape[0]):
    ft[i] = int(final[initial_to_final[i]][8]) == 0 and lib[i] == 0
ax[0].scatter(initial[ft, 0], initial[ft, 1], s=1, c="k", label="survived but not librating")
for i in range(final.shape[0]):
    ft[i] = int(final[i][8]) == 0 and lib[final_to_initial[i]] == 0
ax[1].scatter(final[ft, 0], final[ft, 1], s=1, c="k", label="survived but not librating")

for i in range(initial.shape[0]):
    ft[i] = int(final[initial_to_final[i]][8]) == 0 and lib[i] == 1
ax[0].scatter(initial[ft, 0], initial[ft, 1], s=1, c="r", label="survived and librating")
for i in range(final.shape[0]):
    ft[i] = int(final[i][8]) == 0 and lib[final_to_initial[i]] == 1
ax[1].scatter(final[ft, 0], final[ft, 1], s=1, c="r", label="survived and librating")


for i in range(initial.shape[0]):
    ft[i] = int(final[initial_to_final[i]][8]) == 2
ax[0].scatter(initial[ft, 0], initial[ft, 1], s=4, c="b", label="dead, outofbounds")
for i in range(final.shape[0]):
    ft[i] = int(final[i][8]) == 2
ax[1].scatter(final[ft, 0], final[ft, 1], s=4, c="b", label="dead, outofbounds")


for i in range(initial.shape[0]):
    val = int(final[initial_to_final[i]][8])
    ft[i] = val == 1153
ax[0].scatter(initial[ft, 0], initial[ft, 1], s=1, c="g", label="dead, neptune")
for i in range(final.shape[0]):
    val = int(final[i][8])
    ft[i] = val == 1153
ax[1].scatter(final[ft, 0], final[ft, 1], s=1, c="g", label="dead, neptune")


for i in range(initial.shape[0]):
    val = int(final[initial_to_final[i]][8])
    ft[i] = val == 897
ax[0].scatter(initial[ft, 0], initial[ft, 1], s=10, c="y", marker="x", label="dead, ur")
for i in range(final.shape[0]):
    val = int(final[i][8])
    ft[i] = val == 897
ax[1].scatter(final[ft, 0], final[ft, 1], s=10, c="y", marker="x", label="dead, ur")

for i in range(initial.shape[0]):
    val = int(final[initial_to_final[i]][8])
    ft[i] = val != 0 and val != 2 and val != 1153 and val != 897
ax[0].scatter(initial[ft, 0], initial[ft, 1], s=20, c="r", marker="x", label="dead, ?")
for i in range(final.shape[0]):
    val = int(final[i][8])
    ft[i] = val != 0 and val != 2 and val != 1153 and val != 897
ax[1].scatter(final[ft, 0], final[ft, 1], s=20, c="r", marker="x", label="dead, ?")


ax[0].plot(np.linspace(50, 70, 1000), 1-30.11/np.linspace(50, 70, 1000))
ax[0].plot(np.linspace(50, 70, 1000), 1-20.5/np.linspace(50, 70, 1000))
ax[1].plot(np.linspace(50, 70, 1000), 1-30.11/np.linspace(50, 70, 1000))
ax[1].plot(np.linspace(50, 70, 1000), 1-20.5/np.linspace(50, 70, 1000))

ax[1].set_title("Final")
ax[0].set_title("Initial")

ax[0].set_xlabel("a (AU)")
ax[1].set_xlabel("a (AU)")
ax[0].set_ylabel("e")
ax[1].set_ylabel("e")






fig, ax = plt.subplots(1, 2)

for i in range(initial.shape[0]):
    ft[i] = int(final[initial_to_final[i]][8]) != 0
ax[0].scatter(initial[ft, 0], initial[ft, 1], s=0.5, c="r", label="non-survivors 1Gyr")

for i in range(initial.shape[0]):
    ft[i] = int(final[initial_to_final[i]][8]) == 0
ax[0].scatter(initial[ft, 0], initial[ft, 1], s=0.5, c="k", label="survivors after 1Gyr")
ax[0].legend()

for i in range(final.shape[0]):
    ft[i] = int(final[i][8]) == 0 and not flib[i] and not flib2[i]
ax[1].scatter(final[ft, 0], final[ft, 1], s=0.5, c="lime", label="non-librators")

for i in range(final.shape[0]):
    ft[i] = int(final[i][8]) == 0 and (flib[i] or flib2[i])
ax[1].scatter(final[ft, 0], final[ft, 1], s=0.5, c="k", label="librators")


X = np.linspace(61.63, 63.63, 300)

ax[1].plot(X, 1-30.11/X, c="k")
ax[0].plot(X, 1-30.11/X, c="k")

ax[1].plot(X, 1-20.7/X, c="k", ls="--")
ax[0].plot(X, 1-20.7/X, c="k", ls="--")


ax[1].legend()
ax[0].legend()

ax[0].set_xlim([61.63,63.63])
ax[1].set_xlim([61.63,63.63])
ax[0].set_ylim([0, 0.7])
ax[1].set_ylim([0, 0.7])

ax[0].set_title("Initial 1Myr librators")
ax[1].set_title("Final 1Gyr survivors")

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

for i in range(initial.shape[0]):
    ft[i] = True
ax[0].scatter(initial[ft, 0], initial[ft, 2] / np.pi * 180, s=1, c="k")

for i in range(final.shape[0]):
    ft[i] = int(final[i][8]) == 0 and (flib[i] or flib2[i])
ax[1].scatter(initial[ft, 1], initial[ft, 2] / np.pi * 180,  s=1, c="k")

ax[0].set_xlim([61.63,63.63])
ax[1].set_xlim([0,0.7])
ax[0].set_ylim([0,50])
ax[1].set_ylim([0,50])

ax[0].set_ylabel("i (deg)")
ax[1].set_ylabel("i (deg)")
ax[0].set_xlabel("a (AU)")
ax[1].set_xlabel("e")

ax[0].tick_params(axis="y",direction="in")
ax[0].tick_params(axis="x",direction="in")
ax[0].set_axisbelow(True)
ax[1].tick_params(axis="y",direction="in")
ax[1].tick_params(axis="x",direction="in")
ax[1].set_axisbelow(True)
ax[0].grid()
ax[1].grid()









plt.show()


A = 62.63
plt.figure()
plt.scatter(initial[lib == 0, 0] - A, initial[lib == 0, 1], c="black", s=1, label="Non-librators")
plt.scatter(initial[lib == 1, 0] - A, initial[lib == 1, 1], c="red", s=1, label="Librators")
plt.scatter(initial[lib2 == 1, 0] - A, initial[lib2 == 1, 1], c="blue", s=1, label="Anti-librators")
plt.xlabel("a")
plt.ylabel("e")
plt.legend()

plt.figure()
#plt.scatter(initial[lib == 0, 0], initial[lib == 0, 2], c="black", s=1, label="Non-librators")
plt.scatter(initial[lib == 1, 0], initial[lib == 1, 2], c="red", s=1, label="Librators")
plt.xlabel("a")
plt.ylabel("i (rad)")

plt.figure()
#plt.scatter(initial[lib == 0, 0], initial[lib == 0, 5], c="black", s=1, label="Non-librators")
plt.scatter(initial[lib == 1, 0], initial[lib == 1, 5], c="red", s=1, label="Librators")
plt.xlabel("a")
plt.ylabel("f initial (rad)")

plt.figure()
#plt.scatter(initial[lib == 0, 1], initial[lib == 0, 2], c="black", s=1, label="Non-librators")
plt.scatter(initial[lib == 1, 1], initial[lib == 1, 2], c="red", s=1, label="Librators")
plt.xlabel("e")
plt.ylabel("i (rad)")

plt.show()
