import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def dense_scatter(ax, x, y, data, ny=400, nx=400, logBar=False, label=None, order=None, defaultValue=np.nan, defaultColor=None, upperLimitColor=None):
    xmax = x.max()
    xmin = x.min()
    ymax = y.max()
    ymin = y.min()

    yedges, xedges = np.linspace(ymin, ymax, ny+1), np.linspace(xmin, xmax, nx+1)

    H = np.empty((nx, ny))
    H[:] = defaultValue
    X, Y = np.meshgrid(xedges, yedges)

    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    if order:
        for i in order:
            xi = x[i]
            yi = y[i]
            dat = data[i]
            xc = int(math.floor((xi - xmin) / dx))
            yc = int(math.floor((yi - ymin) / dy))

            if yc == ny: yc -= 1
            if xc == nx: xc -= 1

            H[xc, yc] = dat
    else:
        for xi, yi, dat in zip(x, y, data):
            xc = int(math.floor((xi - xmin) / dx))
            yc = int(math.floor((yi - ymin) / dy))

            if yc == ny: yc -= 1
            if xc == nx: xc -= 1

            H[xc, yc] = dat


    H = H.T[::-1, :]
    ax.grid(False)

    vmax = data.max()

    cm = matplotlib.cm.get_cmap("viridis")

    if defaultColor:
        cm.set_bad(color=defaultColor)

    if np.isnan(defaultValue):
        H = np.ma.masked_array(H, np.isnan(H))

    if upperLimitColor:
        print(vmax)
        vmax = np.nextafter(vmax, np.NINF)
        vmax -= 1
        print(vmax)
        cm.set_over(color=upperLimitColor)

    if logBar:
        norm = matplotlib.colors.LogNorm(vmax=vmax)
    else:
        norm = matplotlib.colors.LinearNorm(vmax=vmax)

    import pdb; pdb.set_trace()

    a = ax.imshow(H, aspect='auto', zorder=0, extent=(xmin, xmax, ymin, ymax), norm=norm, cmap=cm)

    cbar = plt.colorbar(a, ax=ax)
    cbar.set_label(label)

def nologHist(ax, x, y, ny=150, nx=150, logBar=True):
    H, xedges, yedges = np.histogram2d(x, y, bins=[nx, ny])

    H = H.T[::-1, :]
    ax.grid(False)
    a = ax.imshow(np.ma.masked_array(H, H == 0), extent=(x.min(), x.max(), y.min(), y.max()), aspect='auto', norm=matplotlib.colors.LogNorm() if logBar else None, cmap="viridis")
    cbar = plt.colorbar(a, ax=ax)
    cbar.set_label("# particles")
    
def logHist(ax, x, y, ymin, ny=150, nx=150, logBar=False):
    diff = np.log(y.max()) - np.log(ymin)
    y_edges = [ymin]
    for i in range(1, ny+1):
        y_edges.append(np.exp(np.log(ymin) + i * diff / ny))

    H, xedges, yedges = np.histogram2d(x, y, bins=[nx, y_edges])
    X, Y = np.meshgrid(xedges, yedges)
    a = ax.pcolormesh(X, Y, np.ma.masked_array(H.T, H.T == 0), norm=matplotlib.colors.LogNorm() if logBar else None)
    cbar = plt.colorbar(a, ax=ax)
    cbar.set_label("# particles")
    

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
