import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

import numpy as np

with open(sys.argv[1]) as f:
	x = []
	y = []
	z = []
	xp = []
	yp = []
	zp = []
	for line in f:
		vals = line.strip().split()
		if (vals[0] == "pl1"):
			xp.append(float(vals[2]))
			yp.append(float(vals[3]))
			zp.append(float(vals[4]))
		else:
			x.append(float(vals[2]))
			y.append(float(vals[3]))
			z.append(float(vals[4]))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z)
ax.plot(xp, yp, zp)
plt.show()
