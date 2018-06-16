#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib
import math
import sys
import numpy as np
import matplotlib.style as style
import numpy as np

def rv2el(m, parts):
	x = parts[0]
	y = parts[1]
	z = parts[2]
	vx = parts[3]
	vy = parts[4]
	vz = parts[5]

	prec = 1e-13
	pi = np.pi

	hx = y*vz - z*vy
	hy = z*vx - x*vz
	hz = x*vy - y*vx
	hsq = hx*hx + hy*hy + hz*hz

	# As long as we are not on a radial orbit, compute elements
	if hsq > prec:
		i = np.arccos(hz / np.sqrt(hsq))

		# compute the longitude of the ascending node */
		if np.abs(i) < prec:
			capom = 0.0
		elif np.abs(pi - np.abs(i)) < prec:
			capom = 0.0
		else:
			capom = np.arctan2(hx, -hy)

		# /* compute some required quantities */
		vsq = vx*vx + vy*vy + vz*vz
		vdotr =  x*vx + y*vy + z*vz
		r = np.sqrt(x*x + y*y + z*z)
		xhat = x/r
		yhat = y/r
		zhat = z/r
		nx = np.cos(capom)
		ny = np.sin(capom)
		mu = m

		# /* compute the Hamilton vector and thus the eccentricity */
		fac = vsq * r - mu
		Px = fac * xhat - vdotr * vx
		Py = fac * yhat - vdotr * vy
		Pz = fac * zhat - vdotr * vz
		modP = np.sqrt( Px*Px + Py*Py + Pz*Pz )
		e = modP / mu

		# /* compute the argument of pericenter */
		if np.abs(e) < prec:
			om = 0.0
		else:
			if (i < prec) or (pi - i < prec):
				om = arctan2(Py,Px)
			else:
				ecosw = (nx*Px + ny*Py)/mu
				om = np.arccos(ecosw / e)
				if np.abs(Pz) > prec:
					# /* resolve sign ambiguity by sign of Pz  */
					om *= np.abs(Pz)/Pz

		# /* compute the orbital energy , and depending on its sign compute
		# the semimajor axis (or pericenter) and true anomaly	  */
		energy = vsq/2.0 - mu/r
		if np.abs(energy) < prec:
			esign = 0		# /* parabolic */
			a = 0.5 * hsq / mu	# /* actually PERICENTRIC DISTANCE */
			if np.abs(vdotr) < prec:
				f = 0.0
			else:
				f = 2.0 * np.arccos(np.sqrt(a/r)) * vdotr/np.abs(vdotr)
		elif energy > 0.0:
			esign = 1 #		/* hyperbolic */
			a = -0.5 * mu/energy  # /* will be negative */

			if np.abs(vdotr) < prec:
				f = 0.0
			else:
				fac =  a * (1.0 - e * e)/r - 1.0
				f =  np.arccos( fac/ e ) * vdotr/np.abs(vdotr)
		else:
			esign = -1 #		/* elliptic */
			a = -0.5 * mu/energy
			if np.abs(e) > prec:
				if np.abs(vdotr) < prec:
					if r < a:		# /* determine apside */
						f = 0.0
					else:
						f = pi
				else:
					fac =  a * (1.0 - e * e)/r - 1.0
					f =  np.arccos( fac/ e ) * vdotr/np.abs(vdotr)
			else:					   # /* compute circular cases */
				fac = (x * nx + y * ny)/r
				f = np.arccos(fac)
				if np.abs(z) > prec:
					# /* resolve sign ambiguity by sign of z  */
					f *= np.abs(z)/z
				elif (i < prec) or (pi - i < prec):
					f = np.arctan2(y,x) * np.cos(i)
	else: 		#		/* PANIC: radial orbit */
		esign = 1			# /* call it hyperbolic */
		a = np.sqrt(x*x + y*y + z*z)
		e = 1/0.
		i = np.asin(z/sqrt(x*x + y*y + z*z) )	#/* latitude above plane */
		capom = np.arctan2(y,x)			#/* azimuth */
		om = 1/0.
		f = 1/0.
	return (a, e, i, capom, om, f)

style.use('ggplot')

smass = 4 * math.pi * math.pi / math.pow(365.25, 2)

with open(sys.argv[2]) as p:
	n = int(p.readline().strip())

	initial = np.zeros((n, 6))
	for i in range(n):
		xyz = list([float(i) for i in p.readline().strip().split()])
		vxyz = list([float(i) for i in p.readline().strip().split()])
		p.readline()
		initial[i, 0] = xyz[0]
		initial[i, 1] = xyz[1]
		initial[i, 2] = xyz[2]
		initial[i, 3] = vxyz[0]
		initial[i, 4] = vxyz[1]
		initial[i, 5] = vxyz[2]

with open(sys.argv[1]) as p:
	n = int(p.readline().strip())

	final = np.zeros((n, 9))
	for i in range(n):
		xyz = list([float(i) for i in p.readline().strip().split()])
		vxyz = list([float(i) for i in p.readline().strip().split()])
		flags = p.readline().strip().split()

		dtime = float(flags[0])
		pid = int(flags[2])

		final[i, 0] = xyz[0]
		final[i, 1] = xyz[1]
		final[i, 2] = xyz[2]
		final[i, 3] = vxyz[0]
		final[i, 4] = vxyz[1]
		final[i, 5] = vxyz[2]
		final[i, 6] = dtime
		final[i, 7] = pid
		final[i, 8] = int(flags[1])
	initial = initial[final[:, 7].astype(np.int32), :]

for i in range(initial.shape[0]):
	initial[i, :6] = rv2el(smass, initial[i, :6])
for i in range(final.shape[0]):
	final[i, :6] = rv2el(smass, final[i, :6])

def plot_orbits():
    a = np.linspace(22, 28, 200)
    plt.plot(a, 1 - 19.2 / a)
    plt.text(26, 0.25, "Perihelion uranus")
    plt.plot(a, 30 / a - 1)
    plt.text(26, 0.15, "Aphelion neptune")


# import pdb; pdb.set_trace()

alive = final[:, 8] > -1
final[final[:, 6] == 0, 6] = final[:, 6].max()

plt.scatter(initial[alive, 0], initial[alive, 1], c=final[alive, 6] / 365.25e6, s=2, norm=matplotlib.colors.LogNorm())
cbar = plt.colorbar()
cbar.set_label('survival time (MYr)')
plt.title('survival time given initial conditions')
plt.xlim([23, 27])
plt.ylim([0.0001, 0.1])
plt.xlabel('a (AU)')
plt.ylabel('e')
plot_orbits()

plt.figure()
plt.scatter(initial[alive, 0], initial[alive, 2] / 2 / np.pi * 360, c=final[alive, 6] / 365.25e6, s=2, norm=matplotlib.colors.LogNorm())
cbar = plt.colorbar()
cbar.set_label('survival time (MYr)')
plt.title('survival time given initial conditions')
plt.xlim([23, 27])
plt.ylim([0, 0.1])
plt.xlabel('a (AU)')
plt.ylabel('i (deg)')

def cFromFlags(flags):
    nflags = []
    for i in range(len(flags)):
        if flags[i] == 0:
            nflags.append("black")
        elif flags[i] == 1:
            nflags.append("red")
        elif flags[i] == 769:
            nflags.append("green")
        elif flags[i] == 1025:
            nflags.append("blue")
        else:
            print(flags[i], "?")
    return nflags

plt.figure()
plt.title('survival time given initial conditions')
plt.scatter(initial[alive, 0], final[alive, 6] / 365.25e6, s=2, c=cFromFlags(final[alive, 8]))
plt.ylim([0.00001, 4500])
plt.xlabel('a (AU)')
plt.yscale('log')
plt.ylabel('survival time (MYr)')

plt.figure()
plt.title('survival time given initial conditions')
plt.scatter(initial[alive, 1], final[alive, 6] / 365.25e6, s=2, c=cFromFlags(final[alive, 8]))
plt.ylim([0.00001, 4500])
plt.xlabel('e')
plt.yscale('log')
plt.ylabel('survival time (MYr)')

plt.figure()
alive = final[:, 8] == 0
plt.title('surviving particle final states')
plt.scatter(final[alive, 0], final[alive, 1], s=5, c="red", label="Initial states")

if len(alive) < 1000:
    for i in range(len(alive)):
        if alive[i]:
            plt.plot([initial[i, 0], final[i, 0]], [initial[i, 1], final[i, 1]], lw=1, c="pink", zorder=-100)
    plt.scatter(initial[alive, 0], initial[alive, 1], s=5, c="blue", label="Final states")

plt.xlabel('a (AU)')
plt.legend()
plt.ylabel('e')
plot_orbits()
plt.show()



