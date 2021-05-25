#!/usr/bin/env python3
import matplotlib
# matplotlib.use("Qt5Agg")

import matplotlib.pyplot as plt
import util
import util2
import matplotlib.colors
import matplotlib
import math
import sys
import numpy as np
import matplotlib.style as style
import numpy as np

# style.use('ggplot')

smass = 4 * math.pi * math.pi / math.pow(365.25, 2)

state, state_int = util2.read_state(sys.argv[1])
id_to_index = {}
for i in range(state_int.shape[0]):
	id_to_index[state_int[i, 0]] = i

data = []
with open(sys.argv[2], "r") as amplfile:
	for line in amplfile:
		spl = line.split()
		if len(spl) == 15:
			ampl = float(spl[14])
			if float(spl[12]) == 1:
				data.append((int(spl[0]), ampl))

a = []
e = []
datanp = np.array(data)
dataind = []

for index, (x, y) in enumerate(data):
	if x not in id_to_index:
		print(x)
		continue
	dataind.append(index)
	a.append(state[id_to_index[x], 0])
	e.append(state[id_to_index[x], 1])

plt.scatter(a, datanp[dataind, 1], c=e, s=0.5)
cbar = plt.colorbar()
cbar.set_label("e")
plt.xlabel("a (AU)")
plt.ylabel("libration amplitude (deg)")

plt.figure()
plt.scatter(a, e, c=datanp[dataind, 1], s=0.5)
cbar = plt.colorbar()
cbar.set_label("libration amplitude")
plt.xlabel("a (AU)")
plt.ylabel("e")
plt.show()
