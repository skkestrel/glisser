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

smass = 4 * math.pi * math.pi / math.pow(365.25, 2)


dc = {}
with open(sys.argv[2]) as p:
    for line in p:
        tokens = [float(x) for x in line.strip().split()]
        dc[int(tokens[0])] = (tokens[1], tokens[2])

with open(sys.argv[1]) as p:
	npl = int(p.readline().strip())
	for i in range(npl * 4): p.readline()

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

util.nologHist(initial[:, 0], emax[initial[:, 7].astype(np.int32)])
plt.xlabel("a (AU)")
plt.ylabel("e_max")

plt.figure()
util.nologHist(initial[:, 0], emax2[initial[:, 7].astype(np.int32)] / 365e6)
plt.xlabel("a (AU)")
plt.ylabel("e_max_t (MYr)")
plt.show()
