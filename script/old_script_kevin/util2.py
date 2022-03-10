import struct
import numpy as np
import math

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
				om = math.atan2(Py,Px)
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
		e = np.inf
		i = math.asin(z/math.sqrt(x*x + y*y + z*z) )	#/* latitude above plane */
		capom = math.atan2(y,x)			#/* azimuth */
		om = np.inf
		f = np.inf
	return (a, e, i, capom, om, f)

def el2rv(mu, parts):
    a = parts[0]
    e = parts[1]
    i = parts[2]
    capom = parts[3]
    om = parts[4]
    f = parts[5]

    prec = 1e-13

    u = om + f
    xhat = math.cos(u)*math.cos(capom) - math.cos(i)*math.sin(capom)*math.sin(u)
    yhat = math.cos(u)*math.sin(capom) + math.cos(i)*math.cos(capom)*math.sin(u)
    zhat = math.sin(i)*math.sin(u)

    hx = math.sin(capom)*math.sin(i)
    hy = -math.cos(capom)*math.sin(i)
    hz = math.cos(i)

    if abs(e - 1) > prec:
        r = a * (1 - e*e) / (1 + e * math.cos(f))
        h = math.sqrt(mu*a*(1.0 - e*e))
    else:
        h = math.sqrt( 2.0*mu*a )
        r = a/( math.cos(f/2.0)*math.cos(f/2.0) )

    x = r * xhat
    y = r * yhat
    z = r * zhat

    thx = hy * zhat - hz * yhat
    thy = hz * xhat - hx * zhat
    thz = hx * yhat - hy * xhat

    thdot =  h/(r*r)
    rdot  =  e*mu*math.sin(f)/h

    vx = r * thdot * thx + rdot * xhat; 
    vy = r * thdot * thy + rdot * yhat; 
    vz = r * thdot * thz + rdot * zhat; 

    return [x, y, z, vx, vy, vz]

def el2rv_M(mu, parts):
    E = ehybrid(parts[1], parts[5])
    f = math.acos((math.cos(E) - parts[1]) / (1 - parts[1] * math.cos(E)))

    E = E - int(E / 2 / math.pi) * 2 * math.pi
    if E > math.pi: E -= 2 * math.pi
    elif E < -math.pi: E += 2 * math.pi

    # fix the sign of f
    if E < 0: f = -f

    return el2rv(mu, [parts[0], parts[1], parts[2], parts[3], parts[4], f])


def get_principal_angle(M):
    M = M.copy()
    for i in range(len(M)):
        M[i] = M[i] - 2 * np.pi * round(M[i] / (2 * np.pi))
        while M[i] < -np.pi:
            M[i] = M[i] + 2 * np.pi
        while M[i] > np.pi:
            M[i] = M[i] - 2 * np.pi
    return M

def ehie(e, m):
    iflag = 0
    print(m)
    nper = int(m / (2 * math.pi))
    m = m - nper * 2 * math.pi
    if m < 0:
        m += 2 * math.pi;
    if m > math.pi:
        m = 2 * math.pi - m;
        iflag = 1
		
    x = math.pow(6 * m, 1. / 3) - m

    for i in range(3):
        sa = math.sin(x + m)
        ca = math.cos(x + m)
        esa = e * sa
        eca = e * ca
        f = x - esa
        fp = 1 - eca
        dx = -f / fp
        dx = -f / (fp + dx * esa / 2)
        dx = -f / (fp + dx * (esa + eca * dx / 3) / 2)
        x = x + dx

    if iflag:
        return 2 * M_PI - m - x
        m = 2 * M_PI - m
    else:
        return m + x

def eget(e, m, n):
    sm = math.sin(m);
    cm = math.cos(m);
    x = m + e * sm * (1 + e * (cm + e * (1 - 1.5 * sm * sm)))

    for i in range(n):
        sx = math.sin(x)
        cx = math.cos(x)
        es = e*sx
        ec = e*cx
        f = x - es - m
        fp = 1 - ec
        fpp = es
        fppp = ec
        dx = -f / fp
        dx = -f / (fp + dx * fpp / 2)
        dx = -f / (fp + dx * fpp / 2 + dx * dx * fppp / 6)
        x = x + dx
    return x
	
def ehybrid(e, m):
    if (e < 0.18): return eget(e, m, 1);
    elif (e < 0.8): return eget(e, m, 2);
    else: return ehie(e, m)

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
                if not np.all(np.isfinite(pl[i, :6])):
                    print("unbound planet {0}".format(i))
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
            if not np.all(np.isfinite(final[i, :6])):
                print("unbound particle {0}".format(i))
            final[i, :6] = rv2el(smass, final[i, :6])

    if sort:
        ind = np.lexsort((final_int[:, 0],))
        final = final[ind, :]
        final_int = final_int[ind, :]

    if read_planets:
        return pl, pl_int, final, final_int
    else:
        return final, final_int


def bsearch(f, npa, partnum, stride):
    base = f.tell()
    left = 0
    right = npa-1

    while left <= right:
        mid = left + (right - left) // 2
        f.seek(base + mid * stride, 0)
        Id, = struct.unpack('<I', f.read(4))

        if Id == partnum:
            ret = f.read(stride - 4)
            f.seek(base, 0)
            return ret
        elif Id > partnum:
            right = mid - 1
        else:
            left = mid + 1

    f.seek(base, 0)
    return None

