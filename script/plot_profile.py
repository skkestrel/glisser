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

style.use('ggplot')

with open('prof/particle-prof.csv') as p:
    X = []
    Y = []
    for line in p:
        tokens = [float(x) for x in line.strip().split(',')]
        npa = tokens[1]
        nstep = tokens[3]
        time = tokens[4]
        Y.append(time / npa / nstep * 100000 * 10000)
        X.append(npa)
    plt.figure()
    plt.yscale('log')
    plt.xscale('log')
    plt.scatter(X, Y)
with open('prof/timeblock-prof.csv') as p:
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
    plt.yscale('log')
    plt.xscale('log')
    plt.scatter(X, Y)
with open('prof/planet-prof.csv') as p:
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
    plt.yscale('log')
    plt.xscale('log')
    plt.scatter(X, Y)
plt.show()
