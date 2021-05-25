#!/usr/bin/env python3

import matplotlib
#matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import math
import sys
import numpy as np
import matplotlib.style as style
import numpy as np
import util

import sys

plt.figure()
with open('swift-profile2.csv') as p:
    X = []
    Y = []
    for line in p:
        tokens = [float(x) for x in line.strip().split(',')]
        # Y.append(time / npa / nstep * 100000 * 10000)
        Y.append(tokens[1])
        X.append(tokens[0])

    plt.scatter(X, Y, s=3, c="g", label="SWIFT (CPU)")
with open('genga-profile.csv') as p:
    X = []
    Y = []
    for line in p:
        tokens = [float(x) for x in line.strip().split(',')]
        # Y.append(time / npa / nstep * 100000 * 10000)
        Y.append(tokens[1])
        X.append(tokens[0])

    plt.scatter(X, Y, s=3, c="r", label="GENGA (K20m)")
"""
with open('v100_prof/particle-prof.csv') as p:
    X = []
    Y = []
    for line in p:
        tokens = [float(x) for x in line.strip().split(',')]
        npa = tokens[1]
        nstep = tokens[3]
        time = tokens[4]
        # Y.append(time / npa / nstep * 100000 * 10000)
        Y.append(time / nstep * 10000)
        X.append(npa)

    plt.scatter(X, Y, s=3, c="y", label="GLISSE (V100)")
"""
with open('1080ti_prof/particle-prof.csv') as p:
    X = []
    Y = []
    for line in p:
        tokens = [float(x) for x in line.strip().split(',')]
        npa = tokens[1]
        nstep = tokens[3]
        time = tokens[4]
        # Y.append(time / npa / nstep * 100000 * 10000)
        Y.append(time / nstep * 10000)
        X.append(npa)

    plt.scatter(X, Y, s=3, c="b", label="GLISSE (1080ti)")
with open('prof-new/particle-prof.csv') as p:
    X = []
    Y = []
    for line in p:
        tokens = [float(x) for x in line.strip().split(',')]
        npa = tokens[1]
        nstep = tokens[3]
        time = tokens[4]
        # Y.append(time / npa / nstep * 100000 * 10000)
        Y.append(time / nstep * 10000)
        X.append(npa)

    plt.scatter(X, Y, s=3, c="k", label="GLISSE (K20m)")
"""
with open('cpu_prof/particle-prof.csv') as p:
    X = []
    Y = []
    for line in p:
        tokens = [float(x) for x in line.strip().split(',')]
        npa = tokens[1]
        nstep = tokens[3]
        time = tokens[4]
        # Y.append(time / npa / nstep * 100000 * 10000)
        Y.append(time / nstep * 10000)
        X.append(npa)

    plt.scatter(X, Y, s=3, c="r", label="GLISSE (CPU)")
"""

plt.grid()

plt.xlabel('Number of particles')
# plt.ylabel('s / 100k particles / 10k timesteps')
plt.ylabel('s / 10k timesteps')
# plt.yscale('log')
# plt.xscale('log')
plt.ylim([0.1, 100])
plt.xlim([-5000, 140000])
plt.yscale("log")
ax = plt.gca()
ax.tick_params(axis="y",direction="in")
ax.tick_params(axis="x",direction="in")
ax.set_axisbelow(True)
plt.legend()
plt.show()
