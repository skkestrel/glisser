import numpy as np
import math
from numpy import linalg as la
import matplotlib.pyplot as plt
from scipy import optimize
import orbitplotkit as opk


def eigen(nmax, mp, apl, ec, inc, Omega, omega, msun=1.0):

    def laplace_integrand(j, s, a, phi):
        return math.cos(j*phi)/((1.-2.*a*math.cos(phi)+a*a)**s)

    def simpson(y0, y1, y2, h):
        a = y0
        b = (3.*(y2-a)-4.*(y2-y1))/h
        c = ((y2-a)/h-b)/h
        return h*(a+h*(0.5*b+h*c/3.))

    def bint(j, s, a):
        NSIMP = 300
        dphi = 2.0*math.pi/float(NSIMP)
        val = 0.0
        for i in range(NSIMP):
            phi0 = i*dphi
            phi1 = phi0+0.5*dphi
            phi2 = phi0+dphi
            y0 = laplace_integrand(j, s, a, phi0)
            y1 = laplace_integrand(j, s, a, phi1)
            y2 = laplace_integrand(j, s, a, phi2)
            val += simpson(y0, y1, y2, dphi)
        return val/math.pi

    def pq_initial_conds(temp):
        j = 0
        t1 = p[0]
        t2 = q[0]
        for i in range(0, nmax):
            t1 = t1-temp[i]*iv[j, i]*np.sin(temp[i+nmax])
            t2 = t2-temp[i]*iv[j, i]*np.cos(temp[i+nmax])
        out = [t1]
        out.append(t2)
        for j in range(1, nmax):
            t1 = p[j]
            t2 = q[j]
            for i in range(0, nmax):
                t1 = t1-temp[i]*iv[j, i]*np.sin(temp[i+nmax])
                t2 = t2-temp[i]*iv[j, i]*np.cos(temp[i+nmax])
            out.append(t1)
            out.append(t2)
        return(out)

    def hk_initial_conds(temp):
        j = 0
        t1 = h[0]
        t2 = k[0]
        for i in range(0, nmax):
            t1 = t1-temp[i]*ev[j, i]*np.sin(temp[i+nmax])
            t2 = t2-temp[i]*ev[j, i]*np.cos(temp[i+nmax])
        out = [t1]
        out.append(t2)
        for j in range(1, nmax):
            t1 = h[j]
            t2 = k[j]
            for i in range(0, nmax):
                t1 = t1-temp[i]*ev[j, i]*np.sin(temp[i+nmax])
                t2 = t2-temp[i]*ev[j, i]*np.cos(temp[i+nmax])
            out.append(t1)
            out.append(t2)
        return(out)

    def inclination(j, t):
        p0 = 0.
        q0 = 0.
        for i in range(nmax):
            p0 += siv[j, i]*np.sin(f[i]*t+gm[i])
            q0 += siv[j, i]*np.cos(f[i]*t+gm[i])
            inc = np.sqrt(p0**2.0+q0**2.0)*180.0/np.pi
        return inc

    def anode(j, t):
        p0 = 0.
        q0 = 0.
        for i in range(nmax):
            p0 += siv[j, i]*np.sin(f[i]*t+gm[i])
            q0 += siv[j, i]*np.cos(f[i]*t+gm[i])
            node = np.arctan2(p0, q0)*180.0/np.pi
        return node

    def eccentricity(j, t):
        h0 = 0.
        k0 = 0.
        for i in range(nmax):
            h0 += sev[j, i]*np.sin(g[i]*t+bt[i])
            k0 += sev[j, i]*np.cos(g[i]*t+bt[i])
            ecc = np.sqrt(h0**2.0+k0**2.0)
        return ecc

    def longperi(j, t):
        h0 = 0.
        k0 = 0.
        for i in range(nmax):
            h0 += sev[j, i]*np.sin(g[i]*t+bt[i])
            k0 += sev[j, i]*np.cos(g[i]*t+bt[i])
            lperi = np.arctan2(h0, k0)*180.0/np.pi
        return lperi
    sev, siv = np.zeros((nmax, nmax)), np.zeros((nmax, nmax))
    ev, iv = np.zeros((nmax, nmax)), np.zeros((nmax, nmax))
    f, g = np.zeros((nmax)), np.zeros((nmax))
    T, gm = np.zeros((nmax)), np.zeros((nmax))
    S, bt = np.zeros((nmax)), np.zeros((nmax))
    p, q, h, k = np.zeros((nmax)), np.zeros(
        (nmax)), np.zeros((nmax)), np.zeros((nmax))

    # input data from a text file that has planet name, planet mass (in solar masses), eccentricity, inclination, longitude of asc. node, and longitude of perihelion (angles in radians)

    lperi = Omega + omega

    p = inc*np.sin(Omega)
    q = inc*np.cos(Omega)
    h = ec*np.sin(lperi)
    k = ec*np.cos(lperi)

    # mean motions
    mm = np.sqrt(4.*np.pi*np.pi*msun/apl**3.)

    Bf = np.zeros((nmax, nmax))
    Af = np.zeros((nmax, nmax))

    alpha = np.zeros((nmax, nmax))
    alphab = np.zeros((nmax, nmax))

    for j in range(0, nmax):
        for i in range(0, nmax):
            if(j != i):
                if(apl[j] < apl[i]):
                    alpha[j, i] = apl[j]/apl[i]
                    alphab[j, i] = apl[j]/apl[i]
                if(apl[i] < apl[j]):
                    alpha[j, i] = apl[i]/apl[j]
                    alphab[j, i] = 1.

    for j in range(0, nmax):
        temp = 0.
        for i in range(0, nmax):
            if(j != i):
                temp = temp+mp[i]/(msun+mp[j])*alpha[j, i] * \
                    alphab[j, i]*bint(1, 1.5, alpha[j, i])
                Af[j, i] = -0.25*mp[i]/(msun+mp[j])*mm[j] * \
                    alpha[j, i]*alphab[j, i]*bint(2, 1.5, alpha[j, i])
                Bf[j, i] = 0.25*mp[i]/(msun+mp[j])*mm[j]*alpha[j, i] * \
                    alphab[j, i]*bint(1, 1.5, alpha[j, i])
        Af[j, j] = mm[j]*0.25*temp
        Bf[j, j] = -mm[j]*0.25*temp

    # unscaled inclination eigenvectors and frequencies
    iv = np.zeros((nmax, nmax))

    # use numpy routines to solve:
    f, iv = la.eig(Bf)

    # unscaled eccentricity eigenvectors and frequencies
    ev = np.zeros((nmax, nmax))
    g = np.zeros((nmax))
    # use numpy routines to solve:
    g, ev = la.eig(Af)

    # print(f*206264.806)
    # print(g*206264.806)

    temp = np.zeros((2*nmax))

    # scale the inclination eigenvectors
    temp0 = np.zeros((2*nmax))
    temp0.fill(0.5)
    for j in range(nmax, 2*nmax):
        temp0[j] = 1.0
    tempout = optimize.broyden2(pq_initial_conds, temp0)
    # print(tempout)
    for j in range(0, nmax):
        T[j] = tempout[j]
        gm[j] = tempout[j+nmax]

    for j in range(nmax):
        for n in range(nmax):
            siv[j, n] = T[n]*iv[j, n]

    temp = np.zeros((2*nmax))

    # scale the eccentricity eigenvectors
    temp1 = np.zeros((2*nmax))
    temp1.fill(0.5)
    for j in range(nmax, 2*nmax):
        temp1[j] = 1.0
    temp1out = optimize.broyden2(hk_initial_conds, temp1)
    for j in range(0, nmax):
        S[j] = temp1out[j]
        bt[j] = temp1out[j+nmax]

    for j in range(nmax):
        for n in range(nmax):
            sev[j, n] = S[n]*ev[j, n]

    return f, g, siv, sev, gm, bt

    # out = open("inclination-frequencies-amplitudes-MD-format.txt","w")
    # out2 = open("eccentricity-frequencies-amplitudes-MD-format.txt","w")

    # for j in range(nmax):
    #     temp = f[j]*180./np.pi*60*60
    #     temp2 = g[j]*180./np.pi*60*60
    #     out.write("%d %f\n" % (j, temp))
    #     out2.write("%d %f\n" % (j, temp2))

    # out.write(" ")
    # out2.write(" ")
    # for j in range(nmax):
    #     out.write("\t %9d " % (j))
    #     out2.write("\t %9d " % (j))
    # out.write("\n")
    # out2.write("\n")

    # for j in range(nmax):
    #     out.write("%d \t" % (j))
    #     out2.write("%d \t" % (j))
    #     for i in range(nmax):
    #         out.write("%10f " % (siv[i,j]*1e5))
    #         out2.write("%10f " % (sev[i,j]*1e5))
    #     out.write("\n")
    #     out2.write("\n")

    # out.close()
    # out2.close()


def plot(nmax):
    tnmax = 2100
    tmax = 1e7  # in years
    dt = tmax/tnmax
    time = np.arange(0, tmax, dt)

    fig = plt.figure(figsize=(12, 8))  # initialize a figure with size (x,y)
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=0.35, hspace=0.35)

    n = 1
    for j in range(0, nmax):
        ax = fig.add_subplot(nmax, 2, n)
        ax.set_ylabel('inclination (deg)')
        ax.set_xlabel('time (yr)')
        inc = [inclination(j, time[i]) for i in range(tnmax)]
        # ax.axis([0,tmax,0,10.])
        ax.scatter(time, inc, marker='o', color='k', s=0.5)
        n = n+2
    n = 2
    for j in range(0, nmax):
        ax2 = fig.add_subplot(nmax, 2, n)
        ax2.set_ylabel('node (deg)')
        ax2.set_xlabel('time (yr)')
        ax2.axis([0, tmax, -180, 180])
        node = [anode(j, time[i]) for i in range(tnmax)]
        ax2.scatter(time, node, marker='o', color='k', s=0.5)
        n = n+2

    plt.savefig('planets-inclination.png')

    fig2 = plt.figure(figsize=(12, 8))  # initialize a figure with size (x,y)
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=0.35, hspace=0.35)

    n = 1
    for j in range(0, nmax):
        ax = fig2.add_subplot(nmax, 2, n)
        ax.set_ylabel('eccentricity')
        ax.set_xlabel('time (yr)')
        ecc = [eccentricity(j, time[i]) for i in range(tnmax)]
        # ax.axis([0,tmax,0,0.2])
        ax.scatter(time, ecc, marker='o', color='k', s=0.5)
        n = n+2
    n = 2
    for j in range(0, nmax):
        ax2 = fig2.add_subplot(nmax, 2, n)
        ax2.set_ylabel('long. of peri (deg)')
        ax2.set_xlabel('time (yr)')
        ax2.axis([0, tmax, -180, 180])
        peri = [longperi(j, time[i]) for i in range(tnmax)]
        ax2.scatter(time, peri, marker='o', color='k', s=0.5)
        n = n+2

    plt.savefig('planets-eccentricity.png')

# np.set_printoptions(precision=15)
# gn = np.array([22.029082923731803,3.688406991256135,2.687768301546138,0.63074111731915])

# data = np.genfromtxt('planets.txt', names=[
#                      'id', 'm', 'a', 'e', 'inc', 'Omega', 'omega', 'ma'])

# ratio = 1.1
# data['a'][2] = data['a'][2] * ratio
# data['a'][3] = data['a'][3] * ratio
# print(data['a'])
# f, g, siv, sev, gamma, beta = eigen(4, data['m'], data['a'], data['e'], data['inc'], data['Omega'], data['omega'])
# g = np.array(g*opk.RAD_TO_ARCSEC)

# print(g)
# g_diff = (g-gn)/gn
# print(g_diff)

