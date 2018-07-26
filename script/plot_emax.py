"""plot_emax.py

Usage:
    plot_emax.py [options] <initial-state> <emax-file>

Options:
    -h, --help                     Show this screen.
    --plot-mmr <pid>               Plot mmr semi-major axes
"""

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
import docopt

args = docopt.docopt(__doc__)

style.use('ggplot')

smass = 4 * math.pi * math.pi / math.pow(365.25, 2)

dc = {}
with open(args["<emax-file>"]) as p:
    p.readline()
    for line in p:
        tokens = [float(x) for x in line.strip().split(',')]
        dc[int(tokens[0])] = (tokens[1], tokens[2])

with open(args["<initial-state>"]) as p:
    npl = int(p.readline().strip())
    pl_a = 0;
    for i in range(npl):
        p.readline()
        xyz = list([float(i) for i in p.readline().strip().split()])
        vxyz = list([float(i) for i in p.readline().strip().split()])
        Id = int(p.readline().strip().split()[0]) 

        if args["--plot-mmr"] and Id == int(args["--plot-mmr"]):
            els = util.rv2el(smass, np.array(xyz + vxyz))
            PLANET_A = els[0]

    n = int(p.readline().strip())

    initial = np.zeros((n, 9))
    for i in range(n):
        xyz = list([float(i) for i in p.readline().strip().split()])
        vxyz = list([float(i) for i in p.readline().strip().split()])
        flags = p.readline().strip().split()

        dtime = float(flags[2])
        pid = int(flags[0])

        initial[i, 0] = xyz[0]
        initial[i, 1] = xyz[1]
        initial[i, 2] = xyz[2]
        initial[i, 3] = vxyz[0]
        initial[i, 4] = vxyz[1]
        initial[i, 5] = vxyz[2]
        initial[i, 6] = dtime
        initial[i, 7] = pid
        initial[i, 8] = int(flags[1])

for i in range(initial.shape[0]):
    initial[i, :6] = util.rv2el(smass, initial[i, :6])

emax = np.zeros(max(dc.keys())+1)
emax2 = np.zeros(max(dc.keys())+1)
for key, value in dc.items():
    emax[key] = value[0]
    emax2[key] = value[1]

mina = initial[:, 0].min()
maxa = initial[:, 0].max()

def f():
    def mmr(a1, deg):
        return (deg * (a1 ** (3/2))) ** (2/3)
    def f2(deg0, deg1):
        a = mmr(PLANET_A, deg0 / deg1)
        if a < maxa and a > mina:
            plt.axvline(x=a)
            plt.text(a, 0, "{0}:{1}".format(deg0, deg1))

    f2(2, 1)
    f2(3, 2)
    f2(4, 3)
    f2(5, 4)

    f2(1, 2)
    f2(2, 3)
    f2(3, 4)
    f2(4, 5)

    f2(3, 5)
    f2(5, 3)

    f2(2, 5)
    f2(5, 2)

util.nologHist(initial[:, 0], emax[initial[:, 7].astype(np.int32)], 150, 300)
plt.xlabel("a_i (AU)")
plt.ylabel("e_max")
if args["--plot-mmr"]:
    f()

plt.figure()
util.nologHist(initial[:, 0], emax2[initial[:, 7].astype(np.int32)] / 365e6, 150, 300)
plt.xlabel("a_i (AU)")
plt.ylabel("e_max_t (MYr)")
if args["--plot-mmr"]:
    f()
plt.show()
