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

mydict = {}

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
        mydict[Id] = i

lib = np.zeros(initial.shape[0])
lib2 = np.zeros(initial.shape[0])

with open(sys.argv[2]) as p:
    p.readline()
    for line in p:
        lib[mydict[int(line.split(',')[0])]] = int(line.split(',')[1])
        lib2[mydict[int(line.split(',')[0])]] = int(line.split(',')[2])

print("rv2el")

for i in range(initial.shape[0]):
    initial[i, :6] = util.rv2el(smass, initial[i, :6])

#plt.scatter(initial[lib == 0, 0], initial[lib == 0, 1], c="black", s=1, label="Non-librators")
plt.scatter(initial[lib == 1, 0], initial[lib == 1, 1], c="red", s=1, label="Librators")
plt.scatter(initial[lib2 == 1, 0], initial[lib2 == 1, 1], c="blue", s=1, label="Anti-librators")
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
