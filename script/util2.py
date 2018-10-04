import numpy as np

def rv2el(m, parts):
	x = parts[0]
	y = parts[1]
	z = parts[2]
	vx = parts[3]
	vy = parts[4]
	vz = parts[5]

	prec = 1e-13
	pi = math.pi

	hx = y*vz - z*vy
	hy = z*vx - x*vz
	hz = x*vy - y*vx
	hsq = hx*hx + hy*hy + hz*hz

	# As long as we are not on a radial orbit, compute elements
	if hsq > prec:
		i = math.acos(hz / math.sqrt(hsq))

		# compute the longitude of the ascending node */
		if abs(i) < prec:
			capom = 0.0
		elif abs(pi - abs(i)) < prec:
			capom = 0.0
		else:
			capom = math.atan2(hx, -hy)

		# /* compute some required quantities */
		vsq = vx*vx + vy*vy + vz*vz
		vdotr =  x*vx + y*vy + z*vz
		r = math.sqrt(x*x + y*y + z*z)
		xhat = x/r
		yhat = y/r
		zhat = z/r
		nx = math.cos(capom)
		ny = math.sin(capom)
		mu = m

		# /* compute the Hamilton vector and thus the eccentricity */
		fac = vsq * r - mu
		Px = fac * xhat - vdotr * vx
		Py = fac * yhat - vdotr * vy
		Pz = fac * zhat - vdotr * vz
		modP = math.sqrt( Px*Px + Py*Py + Pz*Pz )
		e = modP / mu

		# /* compute the argument of pericenter */
		if abs(e) < prec:
			om = 0.0
		else:
			if (i < prec) or (pi - i < prec):
				om = atan2(Py,Px)
			else:
				ecosw = (nx*Px + ny*Py)/mu
				om = math.acos(ecosw / e)
				if abs(Pz) > prec:
					# /* resolve sign ambiguity by sign of Pz  */
					om *= abs(Pz)/Pz

		# /* compute the orbital energy , and depending on its sign compute
		# the semimajor axis (or pericenter) and true anomaly	  */
		energy = vsq/2.0 - mu/r
		if abs(energy) < prec:
			esign = 0		# /* parabolic */
			a = 0.5 * hsq / mu	# /* actually PERICENTRIC DISTANCE */
			if abs(vdotr) < prec:
				f = 0.0
			else:
				f = 2.0 * math.acos(math.sqrt(a/r)) * vdotr/abs(vdotr)
		elif energy > 0.0:
			esign = 1 #		/* hyperbolic */
			a = -0.5 * mu/energy  # /* will be negative */

			if abs(vdotr) < prec:
				f = 0.0
			else:
				fac =  a * (1.0 - e * e)/r - 1.0
				f =  math.acos( fac/ e ) * vdotr/abs(vdotr)
		else:
			esign = -1 #		/* elliptic */
			a = -0.5 * mu/energy
			if abs(e) > prec:
				if abs(vdotr) < prec:
					if r < a:		# /* determine apside */
						f = 0.0
					else:
						f = pi
				else:
					fac =  a * (1.0 - e * e)/r - 1.0
					f =  math.acos( fac/ e ) * vdotr/abs(vdotr)
			else:					   # /* compute circular cases */
				fac = (x * nx + y * ny)/r
				f = math.acos(fac)
				if abs(z) > prec:
					# /* resolve sign ambiguity by sign of z  */
					f *= abs(z)/z
				elif (i < prec) or (pi - i < prec):
					f = math.atan2(y,x) * math.cos(i)
	else: 		#		/* PANIC: radial orbit */
		esign = 1			# /* call it hyperbolic */
		a = math.sqrt(x*x + y*y + z*z)
		e = 1/0.
		i = math.asin(z/sqrt(x*x + y*y + z*z) )	#/* latitude above plane */
		capom = math.atan2(y,x)			#/* azimuth */
		om = 1/0.
		f = 1/0.
	return (a, e, i, capom, om, f)


def get_principal_angle(M):
    M = M.copy()
    for i in range(len(M)):
        M[i] = M[i] - 2 * np.pi * round(M[i] / (2 * np.pi))
        while M[i] < -np.pi:
            M[i] = M[i] + 2 * np.pi
        while M[i] > np.pi:
            M[i] = M[i] - 2 * np.pi
    return M

def get_M(data):
    E = np.arccos((data[1] + np.cos(data[5])) / (1 + data[1] * np.cos(data[5])))
    E = E * np.sign(data[5])
    M = E - data[1] * np.sin(E)
    return M

def read_state(filename, to_elements=True, read_planets=False, sort=True):
    with open(filename) as p:
        npl = int(p.readline().strip())
        pl = np.zeros((npl, 7))
        pl_int = np.zeros((npl, 1), dtype=np.int32)
        pl_a = 0;

        for i in range(npl):
            m = float(p.readline().strip())
            if i == 0:
                smass = m

            xyz = list([float(i) for i in p.readline().strip().split()])
            vxyz = list([float(i) for i in p.readline().strip().split()])

            Id = int(p.readline().strip().split()[0]) 
            pl_int[i, 0] = Id
            if xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2] < 1e-14:
                pl[i, :6] = 0
            else:
                pl[i, :6] = rv2el(smass, np.array(xyz + vxyz))
            pl[i, 6] = m

        n = int(p.readline().strip())
        final = np.zeros((n, 7))
        final_int = np.zeros((n, 2), dtype=np.int32)
        for i in range(n):
            xyz = list([float(i) for i in p.readline().split()])
            vxyz = list([float(i) for i in p.readline().split()])
            flags = p.readline().strip().split()

            dtime = float(flags[2])
            pid = int(flags[0])

            final[i, 0] = xyz[0]
            final[i, 1] = xyz[1]
            final[i, 2] = xyz[2]
            final[i, 3] = vxyz[0]
            final[i, 4] = vxyz[1]
            final[i, 5] = vxyz[2]
            final[i, 6] = dtime
            final_int[i, 0] = pid
            final_int[i, 1] = int(flags[1])

    if to_elements:
        print("rv2el") 
        for i in range(final.shape[0]):
            final[i, :6] = rv2el(smass, final[i, :6])

    if sort:
        ind = np.lexsort((final_int[:, 0],))
        final = final[ind, :]
        final_int = final_int[ind, :]

    if read_planets:
        return pl, pl_int, final, final_int
    else:
        return final, final_int

