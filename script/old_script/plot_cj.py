"""plot_cj.py

Usage:
    plot_cj.py [options] <stat> <cjmet>

Options:
    -h, --help                     Show this screen.
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

# style.use('ggplot')

stat = {}

with open(args["<stat>"]) as p:
    content = p.readlines()
    for i in range(0, len(content), 3):
        ind = int(content[i].strip())
        mind = float(content[i + 2].strip().split(' ')[4])
        stat[ind] = mind

cj = {}

with open(args["<cjmet>"]) as p:
    p.readline()
    for line in p:
        tokens = [float(x) for x in line.strip().split(',')]
        cj[int(tokens[0])] = (tokens[1], tokens[2], tokens[3])

ids = []
mind = []
cj1 = []
cj2 = []
cj3 = []

for key, val in stat.items():
    ids.append(key)
    mind.append(val)
    cj1.append(cj[key][0])
    cj2.append(cj[key][1])
    cj3.append(cj[key][2])



plt.figure()
plt.scatter(mind, cj2)
plt.xlabel("min dist")
plt.ylabel("max cj err")


plt.figure()
plt.scatter(cj3, cj2)
plt.xlabel("da")
plt.ylabel("max cj err")


plt.figure()
plt.scatter(cj3, cj1)
plt.xlabel("da")
plt.ylabel("final cj err")

plt.figure()
plt.scatter(cj2, cj1)
plt.xlabel("max err")
plt.ylabel("final cj err")

plt.show()
