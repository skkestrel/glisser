"""plot_emax.py

Usage:
    plot_emax.py [options] <initial-state> <emax-file>

Options:
    -h, --help                     Show this screen.
    --plot-mmr <pid>               Plot mmr semi-major axes (can specify list of planets)
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

planet_ids = list([int(x) for x in args['--plot-mmr'].split(',')]) if args['--plot-mmr'] else []
planet_as = {}

with open(args["<initial-state>"]) as p:
    npl = int(p.readline().strip())
    pl_a = 0;
    for i in range(npl):
        p.readline()
        xyz = list([float(i) for i in p.readline().strip().split()])
        vxyz = list([float(i) for i in p.readline().strip().split()])
        Id = int(p.readline().strip().split()[0]) 

        if Id in planet_ids:
            els = util.rv2el(smass, np.array(xyz + vxyz))
            planet_as[Id] = els[0]

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

styles = ['-', '--', '-.', ':']

def f(ax):
    def mmr(a1, deg):
        return (deg * (a1 ** (3/2))) ** (2/3)
    def f2(deg0, deg1):
        for index, i in enumerate(planet_ids):
            a = mmr(planet_as[i], deg0 / deg1)
            if a < maxa and a > mina:
                trans = matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)

                ax.axvline(x=a, ls=styles[index % len(styles)])
                ax.text(a, 1 + 0.05 * index, "{0}:{1}".format(deg0, deg1), transform=trans)

    for i in range(1, 50):
        f2(i, i+1)
        f2(i+1, i)
    for i in range(1, 50):
        if i % 2 == 0: continue
        f2(i, i+2)
        f2(i+2, i)
    for i in range(1, 50):
        if i % 3 == 0: continue
        f2(i, i+3)
        f2(i+3, i)
    for i in range(1, 50):
        if i % 2 == 0: continue
        f2(i, i+4)
        f2(i+4, i)
        '''
    for i in range(1, 50):
        if i % 5 == 0: continue
        f2(i, i+5)
        f2(i+5, i)
    for i in range(1, 20):
        if i % 2 == 0: continue
        if i % 3 == 0: continue
        f2(i, i+6)
        f2(i+6, i)
    for i in range(1, 20):
        if i % 7 == 0: continue
        f2(i, i+7)
        f2(i+7, i)
    for i in range(1, 20):
        if i % 2 == 0: continue
        f2(i, i+8)
        f2(i+8, i)
        '''

fig, ax = plt.subplots(2, 2)

util.nologHist(ax[0, 0], initial[:, 0], emax[initial[:, 7].astype(np.int32)], 150, 300, False)
ax[0, 0].set_xlabel("a_i (AU)")
ax[0, 0].set_ylabel("e_max")
f(ax[0, 0])

util.nologHist(ax[1, 0], initial[:, 0], emax2[initial[:, 7].astype(np.int32)] / 365e6, 150, 300, False)
ax[1, 0].set_xlabel("a_i (AU)")
ax[1, 0].set_ylabel("e_max_t (MYr)")
f(ax[1, 0])

util.dense_scatter(ax[0, 1], initial[:, 0], emax[initial[:, 7].astype(np.int32)], initial[:, 1], label="e_i")
ax[0, 1].set_xlabel("a_i (AU)")
ax[0, 1].set_ylabel("e_max")
f(ax[0, 1])

if len(planet_as) > 0:
    for index, Id in enumerate(planet_ids):
        ax[0, 1].set_prop_cycle(None)
        ax[0, 1].plot([], [], label="Planet {0}".format(Id), ls=styles[index % len(styles)], transform=ax[0, 1].transAxes)
    ax[0, 1].legend()

util.dense_scatter(ax[1, 1], initial[:, 0], emax2[initial[:, 7].astype(np.int32)] / 365e6, initial[:, 1], label="e_i")
ax[1, 1].set_xlabel("a_i (AU)")
ax[1, 1].set_ylabel("e_max_t (MYr)")
f(ax[1, 1])

plt.show()
