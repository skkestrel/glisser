import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import orbitplotkit as opk
import rebound
from scipy.special import ellipk
import time


NUM = 360
EPS = 1e-7


def disturbFunction(rj, r):
    '''
    rj - position vector of the planet
    r - position vector of the small object
    '''
    # rj_dist = np.linalg.norm(rj)
    delta_dist = np.linalg.norm(rj-r)
    # print(rj_dist, delta_dist)

    result = 1/delta_dist
    return -result

def disturbFunctionRing(aj, r_dist, x2y2):
    '''
    aj - semi-major axis of the planet
    r_dist - distance from the center
    x2y2 - x^2 + y^2
    '''

    temp1 = aj*np.sqrt(x2y2)
    temp2 = r_dist**2 + aj**2

    m = 4*temp1/(temp2 + 2*temp1)
    result = 4/np.sqrt(temp2) * np.sqrt(1-m/2) * ellipk(m)
    return result

def doubleAverage(mj, aj, a, G, H, g, num):
    if G > 0:
        cosi = H/G
        if cosi >= 1:
            inc = 0
        else:
            inc = np.arccos(cosi)
    else:
        inc = 0
    e = np.sqrt(1 - G**2/a)

    sim = rebound.Simulation()
    # sim.G = 4*np.pi**2
    sim.add(m=1.)
    sim.add(m=mj, a=aj, e=1e-7, inc=1e-7)
    sim.add(m=0, a=a, e=e, inc=inc, Omega=0, omega=g)

    sun = sim.particles[0]
    planet = sim.particles[1]
    particle = sim.particles[2]

    particle.M=0
    orbit = particle.calculate_orbit(sun)

    num_j = num
    num_i = num
    num_total = num_i * num_j
    rj_list = []
    r_list = []
    for j in np.arange(num_j):
        Mj = np.deg2rad(j/num_j * 360)
        planet.omega = 0
        planet.Omega = 0
        planet.M = Mj
        orbit_j = planet.calculate_orbit(sun)
        rj = np.array([planet.x, planet.y, planet.z])
        rj_list.append(rj)
    for i in np.arange(num_i):
        M = np.deg2rad(i/num_i * 360)
        particle.omega = g
        particle.Omega = 0
        particle.M = M
        orbit = particle.calculate_orbit(sun)
        r = np.array([particle.x, particle.y, particle.z])
        r_list.append(r)

    result = 0
    count = 0
    for rj in rj_list:
        for r in r_list:
            result += disturbFunction(rj, r)
            count += 1
    return result/num_total*mj


def doubleAverageRing(mj, aj, a, G, H, g, num):
    if G > 0:
        cosi = H/G
        # inc = np.arccos(cosi)
        sini = np.sqrt(1-cosi**2)
    else:
        inc = 0
    e = np.sqrt(1 - G**2/a)

    if a < 0 or e < 0 or e > 1:
        return np.nan



    # sim = rebound.Simulation()
    # sim.G = 4*np.pi**2
    # sim.add(m=1.)
    # sim.add(m=mj, a=aj, e=1e-7, inc=1e-7)
    # sim.add(m=0, a=a, e=e, inc=inc, Omega=0, omega=g)

    # sun = sim.particles[0]
    # planet = sim.particles[1]
    # particle = sim.particles[1]

    # particle.M=0
    # orbit = particle.calculate_orbit(sun)
    result = 0
    # count = 0
    for i in np.arange(num):
        M = np.deg2rad(i/num * 360)
        f = opk.M2F(M, e)
        r_dist = a*(1-e**2)/(1+e*np.cos(f))
        x2y2 = r_dist**2 * (1 - np.sin(g+f)**2 * sini**2)


        # particle.omega = g
        # particle.Omega = 0
        # particle.M = M
        # orbit = particle.calculate_orbit(sun)
        # r = np.array([particle.x, particle.y, particle.z])
        result += disturbFunctionRing(aj, r_dist, x2y2)
        # count += 1
    return result/num*mj

def dhdt(mj, aj, a, e, inc, g, num = NUM, epsilon = EPS):
    L = np.sqrt(a)
    G = L*np.sqrt(1-e**2)
    H = G*np.cos(inc)

    Y2 = doubleAverage(mj, aj, a, G, H+epsilon, g, num)
    Y1 = doubleAverage(mj, aj, a, G, H-epsilon, g, num)

    result = (Y2-Y1)/(2*epsilon)
    return -2*np.pi*result

def dhdtRing(mj, aj, a, e, inc, g, num = NUM, epsilon = EPS):
    L = np.sqrt(a)
    G = L*np.sqrt(1-e**2)
    H = G*np.cos(inc)

    Y2 = doubleAverageRing(mj, aj, a, G, H+epsilon, g, num)
    Y1 = doubleAverageRing(mj, aj, a, G, H-epsilon, g, num)

    result = (Y2-Y1)/(2*epsilon)
    return result

def dgdtRing(mj, aj, a, e, inc, g ,num = NUM, epsilon = EPS):
    L = np.sqrt(a)
    G = L*np.sqrt(1-e**2)
    H = G*np.cos(inc)

    Y2 = doubleAverageRing(mj, aj, a, G+epsilon, H, g, num)
    Y1 = doubleAverageRing(mj, aj, a, G-epsilon, H, g, num)

    result = (Y2-Y1)/(2*epsilon)
    return result

def dpidtRing(mj, aj, a, e, inc, g ,num = NUM, epsilon = EPS):
    L = np.sqrt(a)
    G = L*np.sqrt(1-e**2)
    H = G*np.cos(inc)

    Y2 = doubleAverageRing(mj, aj, a, G, H+epsilon, g, num)
    Y1 = doubleAverageRing(mj, aj, a, G, H-epsilon, g, num)

    result = (Y2-Y1)/(2*epsilon)

    Y2 = doubleAverageRing(mj, aj, a, G+epsilon, H, g, num)
    Y1 = doubleAverageRing(mj, aj, a, G-epsilon, H, g, num)

    result += (Y2-Y1)/(2*epsilon)

    return result

# end = time.time()
# print(end - start)

# start = time.time()
# for i in np.arange(3):
#     r2 = dhdtRing(mj, aj, a, e, inc, g, NUM)
# end = time.time()
# print(end - start)

# print(r1, r2 ,r2/r1)

# a = 40
# e = 0.2
# inc = np.deg2rad(10)

# print(a,e,inc)
# # print(L,G,H)

# print(dhdtRing(1e-3, 30, a, e, inc, 0, epsilon=EPS))
# print(dpidtRing(1e-3, 30, a, e, inc, 0, epsilon=EPS))
# x = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
# y = []
# for eps in x:
#     result = )
#     print(eps, result)
#     y.append(result+2.3375567375555405e-06)
# print(y)

# fig, ax = plt.subplots()
# ax.plot(x,y)
# ax.set_xscale('log')

# fig.savefig("test_out.png", dpi=300)