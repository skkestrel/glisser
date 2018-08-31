#!/usr/bin/env python3

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import math
import sys
import numpy as np
import matplotlib.style as style
import numpy as np
import util

import sys

# style.use('ggplot')

Dir = sys.argv[1]

with open('{0}/particle-prof.csv'.format(Dir)) as p:
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
    plt.figure()
    plt.grid()

    plt.xlabel('Number of particles')
    # plt.ylabel('s / 100k particles / 10k timesteps')
    plt.ylabel('s / 10k timesteps')
    # plt.yscale('log')
    # plt.xscale('log')
    plt.ylim([-0.05, 3])
    plt.xlim([-1000, 140000])
    plt.scatter(X, Y, s=3, c="k")
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
with open('{0}/timeblock-prof.csv'.format(Dir)) as p:
    X = []
    Y = []
    for line in p:
        tokens = [float(x) for x in line.strip().split(',')]
        npa = tokens[1]
        nblock = tokens[2]
        nstep = tokens[3]
        time = tokens[4]
        Y.append(time / npa / nstep * 100000 * 10000)
        X.append(nblock)
    plt.figure()
    plt.grid()
    plt.xlabel('Timeblock size')
    plt.ylabel('s / 100k particles / 10k timesteps')
    plt.yscale('log')
    plt.xscale('log')
    plt.scatter(X, Y, s=3, c="k")
    ax = plt.gca()
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.set_axisbelow(True)
with open('{0}/planet-prof.csv'.format(Dir)) as p:
    X = []
    Y = []
    for line in p:
        tokens = [float(x) for x in line.strip().split(',')]
        npl = tokens[0]
        npa = tokens[1]
        nstep = tokens[3]
        time = tokens[4]
        Y.append(time / npa / nstep * 100000 * 10000)
        X.append(npl)
    plt.figure()
    plt.xlabel('Number of planets')
    plt.ylabel('s / 10k timesteps')
    #plt.yscale('log')
    #plt.xscale('log')
    plt.scatter(X, Y)
plt.show()
