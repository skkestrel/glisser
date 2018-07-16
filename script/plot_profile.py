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

with open('particle-prof.csv') as p:
	X = []
	Y = []
	for line in p:
		tokens = [float(x) for x in line.strip().split(',')]
		npa = tokens[1]
		nstep = tokens[3]
		time = tokens[4]
		Y.append(time)
		X.append(npa / nstep * 100000)
	plt.plot(X, Y)
	plt.show()
