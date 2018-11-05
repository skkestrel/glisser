"""plot_emax.py

Usage:
    plot_emax.py [options] <initial-state> <final-state> <emax-file>

Options:
    -h, --help                     Show this screen.
    --plot-mmr <pid>               Plot mmr semi-major axes (can specify list of planets)
"""

import matplotlib
# matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import util2
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

pl, pl_int, initial, initial_int = util2.read_state(args["<initial-state>"], read_planets=True)
final, final_int = util2.read_state(args["<final-state>"])

if not np.equal(final_int[:, 0], initial_int[:, 0]).all():
    print("?")
    sys.exit(-1)

for i in range(pl.shape[0]):
    if pl_int[i, 0] in planet_ids:
        planet_as[pl_int[i, 0]] = pl[i, 0]
#planet_as[4] = 30.1098
#planet_as[3] = 19.2183
planet_as[4] = 29.9871
planet_as[3] = 19.3173
print(planet_as)

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
final[final[:, 6] == 0, 6] = final[:, 6].max() * 1.1

util.nologHist(ax[0, 0], initial[:, 0], emax[initial_int[:, 0]], 150, 300, False)
ax[0, 0].set_xlabel("a_i (AU)")
ax[0, 0].set_ylabel("e_max")
f(ax[0, 0])

# util.nologHist(ax[1, 0], initial[:, 0], emax2[initial_int[:, 0]] / 365e6, 150, 300, False)
util.dense_scatter(ax[1, 0], initial[:, 0], emax2[initial_int[:, 0]] / 365e6, emax[initial_int[:, 0]], label="e_max", ny=400, nx=400, vmax=0.32)
ax[1, 0].set_xlabel("a_i (AU)")
ax[1, 0].set_ylabel("e_max_t (MYr)")
f(ax[1, 0])

util.dense_scatter(ax[0, 1], initial[:, 0], emax[initial_int[:, 0]], initial[:, 1], label="e_i")
ax[0, 1].set_xlabel("a_i (AU)")
ax[0, 1].set_ylabel("e_max")
f(ax[0, 1])

if len(planet_as) > 0:
    for index, Id in enumerate(planet_ids):
        ax[0, 1].set_prop_cycle(None)
        ax[0, 1].plot([], [], label="Planet {0}".format(Id), ls=styles[index % len(styles)], transform=ax[0, 1].transAxes)
    ax[0, 1].legend()

# util.dense_scatter(ax[1, 1], initial[:, 0], emax2[initial_int[:, 0]] / 365e6, initial[:, 1], label="e_i")
util.dense_scatter(ax[1, 1], initial[:, 0], emax2[initial_int[:, 0]] / 365e6, final[:, 6] / 365e6, label="deathtime (MYr)", upperLimitColor="r")
ax[1, 1].set_xlabel("a_i (AU)")
ax[1, 1].set_ylabel("e_max_t (MYr)")
f(ax[1, 1])

fig, ax = plt.subplots(2, 1)
final[final[:, 6] == 0, 6] = final[:, 6].max() * 1.1

# util.nologHist(ax[1, 0], initial[:, 0], emax2[initial_int[:, 0]] / 365e6, 150, 300, False)
util.dense_scatter(ax[0], initial[:, 0], emax[initial_int[:, 0]], initial[:, 1], label="e_i")
ax[0].set_xlabel("a_i (AU)")
ax[0].set_ylabel("e_max")

# util.dense_scatter(ax[1, 1], initial[:, 0], emax2[initial_int[:, 0]] / 365e6, initial[:, 1], label="e_i")
util.dense_scatter(ax[1], initial[:, 0], emax2[initial_int[:, 0]] / 365e6, final[:, 6] / 365e6, label="deathtime (MYr)", upperLimitColor="crimson")
ax[1].set_xlabel("a_i (AU)")
ax[1].set_ylabel("e_max_t (MYr)")

plt.show()
