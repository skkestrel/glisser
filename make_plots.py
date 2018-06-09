import matplotlib.pyplot as plt
import numpy as np
import math

with open("ics.out") as p:
	n = int(p.readline().strip())
	a_s = []
	t_s = []
	for i in range(n):
		xyz = list([float(i) for i in p.readline().split()])
		vxyz = list([float(i) for i in p.readline().split()])
		r = math.sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2])
		f = p.readline()
		t_s.append(float(f.split()[0]))
		a_s.append(r)
	a_s = np.array(a_s)
	t_s = np.array(t_s)

	t_s[t_s == 0] = t_s.max()
	plt.scatter(a_s, t_s)
	plt.show()
