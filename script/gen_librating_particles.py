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

pl, pl_int, initial, initial_int = util2.read_state(sys.argv[1], read_planets=True)

lib = set()

with open(sys.argv[2]) as p:
    p.readline()
    for line in p:
        if int(line.split(',')[1]) == 1 or int(line.split(',')[2]) == 1:
             lib.add(int(line.split(',')[0]))

with open("5-2librators.txt", "w") as f:
	f.write("-5 {0} {1} {2} {3} {4} {5}\n".format(*pl[4, :5], util2.get_M(pl[4, :6])))
	for i in range(initial.shape[0]):
		if initial_int[i, 0] in lib:
			f.write("{6} {0} {1} {2} {3} {4} {5}\n".format(*initial[i, :5], util2.get_M(initial[i, :6]), initial_int[i, 0]))
