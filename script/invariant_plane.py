import numpy as np
import matplotlib.pyplot as plt
import orbitplotkit as opk
from pylaplace import LaplaceCoefficient
import linear_secular_solution as lss
import double_average as da
from scipy.optimize import root_scalar
from scipy.optimize import brenth
import time


plt.style.use('aiur')

def alpha_j(a, aj):
    if aj <= a:
        return aj/a
    else:
        return a/aj


def alpha_jb(a, aj):
    if aj <= a:
        return 1
    else:
        return a/aj


def meanMotion(a):
    """
    unit: rad/year
    """
    mu = 4*np.pi**2
    return np.sqrt(mu/a**3)

def Bj_func(a, aj, mj, msun):
    laplace = LaplaceCoefficient(method='Brute')
    b32 = laplace(alpha_j(a, aj), 1.5, 1, 1, 1)
    result = meanMotion(a)/4 * mj / msun * alpha_j(a, aj) * \
            alpha_jb(a, aj) * b32
    return result


class InvariantPlane(object):

    msun = 1.0

    def __init__(self, nmax, mp, apl, ec, inc, Omega, omega):
        self.nmax = nmax
        self.mp = mp
        self.apl = apl

        self.f, self.g, self.siv, self.sev, self.gamma, self.beta = lss.eigen(
            nmax, mp, apl, ec, inc, Omega, omega, msun=1.0)

    def printdata(self):
        print(self.mp)
        print(self.apl)
        print(self.f)
        print(self.siv)
        print(self.gamma)
        print(self.nmax)

    def Bj_func(self, a, aj, mj):
        laplace = LaplaceCoefficient(method='Brute')
        b32 = laplace(alpha_j(a, aj), 1.5, 1, 1, 1)
        result = meanMotion(a)/4 * mj / self.msun * alpha_j(a, aj) * \
            alpha_jb(a, aj) * b32
        return result

    def Bj_func_da(self, a, aj, mj, e, inc, g):
        if aj < 25:
            result = da.dhdtRing(mj, aj, a, e, inc, g, num = 10)
        else:
            result = da.dhdtRing(mj, aj, a, e, inc, g, num = 180)
        return result

    def B_func(self, a):
        result = 0
        for aj, mj in zip(self.apl, self.mp):
            result -= self.Bj_func(a, aj, mj)
        return result

    def B_func_da(self, a, e, inc, g):
        result = 0
        for aj, mj in zip(self.apl, self.mp):
            result -= self.Bj_func_da(a, aj, mj, e, inc, g)
        return result

    def Aj_func_da(self, a, aj, mj, e, inc, g):
        if aj < 25:
            result = da.dpidtRing(mj, aj, a, e, inc, g, num = 10)
        else:
            result = da.dpidtRing(mj, aj, a, e, inc, g, num = 180)
        return result


    def A_func_da(self, a, e, inc, g):
        result = 0
        for aj, mj in zip(self.apl, self.mp):
            result -= self.Aj_func_da(a, aj, mj, e, inc, g)
        return result

    def mu_func(self, a, sivi):
        result = 0
        for aj, mj, Iji in zip(self.apl, self.mp, sivi):
            result += self.Bj_func(a, aj, mj) * Iji
        return result

    def mu_func_da(self, sivi, Bjs):
        result = 0
        for Iji, Bj in zip(sivi, Bjs):
            result += Bj * Iji
        return result

    def q_p_forced(self, t, a):
        q0 = 0
        p0 = 0
        Bfunc = self.B_func(a)
        for i, fi, gammai in zip(np.arange(self.nmax), self.f, self.gamma):
            sivi = self.siv[:, i]
            q0 -= self.mu_func(a, sivi)/(Bfunc - fi) * \
                np.cos(fi*t + gammai)
            p0 -= self.mu_func(a, sivi) / \
                (Bfunc - fi) * np.sin(fi*t + gammai)
        # print(Bfunc*opk.RAD_TO_ARCSEC)
        return q0, p0

    def q_p_forced_da(self, t, a, e, inc, g):
        q0 = 0
        p0 = 0
        Bj = []
        Bfunc = 0
        for aj, mj in zip(self.apl, self.mp):
            temp = self.Bj_func_da(a, aj, mj, e, inc, g)
            Bj.append(temp)
            Bfunc -= temp
        # print(np.array(Bj)*opk.RAD_TO_ARCSEC)
        # print(Bfunc*opk.RAD_TO_ARCSEC)
        # print(np.array(self.f)*opk.RAD_TO_ARCSEC)
        for i, fi, gammai in zip(np.arange(self.nmax), self.f, self.gamma):
            sivi = self.siv[:, i]
            mufunc = self.mu_func_da(sivi, Bj)
            q0 -= mufunc/(Bfunc - fi) * np.cos(fi*t + gammai)
            p0 -= mufunc/(Bfunc - fi) * np.sin(fi*t + gammai)
        return q0, p0

    def q_p_forced_freq(self, t, a, freq):
        """
        input: freq (rad per year)
        """
        q0 = 0
        p0 = 0
        Bfunc = freq
        for i, fi, gammai in zip(np.arange(self.nmax), self.f, self.gamma):
            sivi = self.siv[:, i]
            q0 -= self.mu_func(a, sivi)/(Bfunc - fi) * \
                np.cos(fi*t + gammai)
            p0 -= self.mu_func(a, sivi) / \
                (Bfunc - fi) * np.sin(fi*t + gammai)
        return q0, p0

    def q_p_full(self, t, a, I, Omega):
        I_free, Omega_free = self.I_Omega_free(a, I, Omega)
        q0, p0 = self.q_p_forced(t, a)
        B = self.B_func(a)
        q = I_free * np.cos(B*t + Omega_free) + q0
        p = I_free * np.sin(B*t + Omega_free) + p0
        return q, p

    def I_Omega_evolution(self, t, a, I, Omega):
        q, p = self.q_p_full(t, a, I, Omega)
        I, Omega = self.I_Omega_trans(q, p)
        return I, Omega

    def q_p_trans(self, I, Omega):
        q = I*np.cos(Omega)
        p = I*np.sin(Omega)
        return q, p

    def I_Omega_trans(self, q, p):
        I = np.sqrt(q**2 + p**2)
        Omega = np.arctan2(p, q)
        return I, Omega

    def I_Omega_forced(self, a):
        """
        input: a (au)
        output: radians
        """
        q0, p0 = self.q_p_forced(0, a)
        I_forced, Omega_forced = self.I_Omega_trans(q0, p0)
        return I_forced, Omega_forced

    def I_Omega_forced_da(self, a, e, inc, g):
        """
        input: a (au)
        output: radians
        """
        q0, p0 = self.q_p_forced_da(0, a, e, inc, g)
        I_forced, Omega_forced = self.I_Omega_trans(q0, p0)
        return I_forced, Omega_forced

    def I_Omega_forced_freq(self, a, freq):
        """
        input: a (au)
        output: radians
        """
        q0, p0 = self.q_p_forced_freq(0, a, freq)
        I_forced, Omega_forced = self.I_Omega_trans(q0, p0)
        return I_forced, Omega_forced

    def I_Omega_free(self, a, I, Omega):
        """
        input: a (au), I (rad), Omega (rad)
        """

        # I_forced, Omega_forced = I_Omega_forced(a)
        # print(np.rad2deg(I_forced), np.rad2deg(Omega_forced))

        q_forced, p_forced = self.q_p_forced(0, a)
        q1, p1 = self.q_p_trans(I, Omega)

        q_free, p_free = q1 - q_forced, p1 - p_forced
        I_free, Omega_free = self.I_Omega_trans(q_free, p_free)

        return I_free, Omega_free

    def I_Omega_free_da(self, a, e, I, Omega, g):
        """
        input: a (au), I (rad), Omega (rad)
        """
        inc = I
        # I_forced, Omega_forced = I_Omega_forced(a)
        # print(np.rad2deg(I_forced), np.rad2deg(Omega_forced))

        q_forced, p_forced = self.q_p_forced_da(0, a, e, inc, g)

        q1, p1 = self.q_p_trans(I, Omega)

        q_free, p_free = q1 - q_forced, p1 - p_forced
        I_free, Omega_free = self.I_Omega_trans(q_free, p_free)

        # print(a,e,np.rad2deg(inc),np.rad2deg(g))
        # print(np.rad2deg(I_free), np.rad2deg(Omega_free))
        # print('')

        return I_free, Omega_free


    def I_Omega_free_qp_forced_da(self, a, e, I, Omega, g):
        """
        input: a (au), e, I (rad), Omega (rad), g (rad), I_ivp (rad), Omega_ivp (rad)
        """

        inc = I
        q_forced, p_forced = self.q_p_forced_da(0, a, e, inc, g)

        q1, p1 = self.q_p_trans(I, Omega)

        q_free, p_free = q1 - q_forced, p1 - p_forced
        I_free, Omega_free = self.I_Omega_trans(q_free, p_free)

        return I_free, Omega_free, q_forced, p_forced

    def I_Omega_free_freq(self, a, I, Omega, freq):
        """
        input: a (au), I (rad), Omega (rad)
        """

        # I_forced, Omega_forced = I_Omega_forced(a)
        # print(np.rad2deg(I_forced), np.rad2deg(Omega_forced))

        q_forced, p_forced = self.q_p_forced_freq(0, a, freq)
        q1, p1 = self.q_p_trans(I, Omega)

        q_free, p_free = q1 - q_forced, p1 - p_forced
        I_free, Omega_free = self.I_Omega_trans(q_free, p_free)

        return I_free, Omega_free


def miniFuncB(a, e, inc):
    data = np.genfromtxt('planets.txt', names=['id', 'mass', 'a', 'e', 'inc', 'Omega', 'omega', 'ma'])
    ivp = InvariantPlane(4, data['mass'], data['a'],
                        data['e'], data['inc'], data['Omega'], data['omega'])
    result = 0
    n = 10
    for g in np.linspace(0,90,n):
        g = np.deg2rad(g)
        result += ivp.B_func_da(a, e, inc, g)
    result /= n
    y = result - opk.f_freqs[7]*opk.ARCSEC_TO_RAD
    return y

def miniFuncB2(e, a, inc, g):
    data = np.genfromtxt('planets.txt', names=['id', 'mass', 'a', 'e', 'inc', 'Omega', 'omega', 'ma'])
    ivp = InvariantPlane(4, data['mass'], data['a'],
                        data['e'], data['inc'], data['Omega'], data['omega'])
    y = ivp.B_func_da(a, e, inc, g) - opk.f_freqs[7]*opk.ARCSEC_TO_RAD
    return y

def findS8(a0, e, inc):
    result= root_scalar(miniFuncB, args=(e, inc), bracket=(38, 44), x0=a0, xtol=1e-10, maxiter=30, method='brentq')
    return result.root, result.converged

def findS8_2(e0, a, inc, g, e1):
    result= root_scalar(miniFuncB2, args=(a, inc, g), bracket=(0.001, e1), x0=e0, xtol=1e-10, maxiter=30, method='brentq')
    return result.root, result.converged


def miniFuncA(a, e, inc):
    data = np.genfromtxt('planets.txt', names=['id', 'mass', 'a', 'e', 'inc', 'Omega', 'omega', 'ma'])
    ivp = InvariantPlane(4, data['mass'], data['a'],
                        data['e'], data['inc'], data['Omega'], data['omega'])
    result = 0
    n = 10
    for g in np.linspace(0,90,n):
        g = np.deg2rad(g)
        result += ivp.A_func_da(a, e, inc, g)
    result /= n
    y = result - opk.g_freqs[7]*opk.ARCSEC_TO_RAD

    return y

def miniFuncA2(e, a, inc):
    data = np.genfromtxt('planets.txt', names=['id', 'mass', 'a', 'e', 'inc', 'Omega', 'omega', 'ma'])
    ivp = InvariantPlane(4, data['mass'], data['a'],
                        data['e'], data['inc'], data['Omega'], data['omega'])
    result = 0
    n = 10
    for g in np.linspace(0,90,n):
        g = np.deg2rad(g)
        result += ivp.A_func_da(a, e, inc, g)
    result /= n
    y = result - opk.g_freqs[7]*opk.ARCSEC_TO_RAD

    return y

def findG8(a0, e, inc):
    result= root_scalar(miniFuncA, args=(e, inc), bracket=(38, 44), x0=a0, xtol=1e-10, maxiter=30, method='brentq')
    return result.root, result.converged

def findG8_2(e0, a, inc, e1):
    result= root_scalar(miniFuncA2, args=(a, inc), bracket=(0.001, e1), x0=e0, xtol=1e-10, maxiter=30, method='brentq')
    return result.root, result.converged

def runINC(e):
    # data = np.genfromtxt('planets.txt', names=[
    #                     'id', 'mass', 'a', 'e', 'inc', 'Omega', 'omega', 'ma'])

    # print(data['a'])

    # ivp = InvariantPlane(4, data['mass'], data['a'],
    #                     data['e'], data['inc'], data['Omega'], data['omega'])
    a0 = 40.4

    # inc_fin = 14
    inc_list = np.logspace(-1,1.5,40)
    # print(inc_list)
    # inc_list[0] = 0.1

    a_list = []
    for inc in inc_list:
        inc_rad = np.deg2rad(inc)
        a, err = findS8(a0, e, inc_rad)
        a_list.append(a)
        print(inc, a, err)
        a0 = a
        if a < 39:
            break


def runECC(inc):
    e0 = 0.18

    # inc_fin = 14
    a_list = np.linspace(39.185, 44, 40)
    q1 = 30.08
    # inc_list[0] = 0.1
    # gs = [45]
    e_list = []
    
    for a in a_list:
        e1 = 1 - q1/a
        e, err  = findG8_2(e0, a, np.deg2rad(inc), e1)
        e_list.append(e)
        e0 = e
        print(a, e, a*(1-e), err)

# runECC(8)
# runINC(0.001)

# data = np.genfromtxt('planets.txt', names=[
#                     'id', 'mass', 'a', 'e', 'inc', 'Omega', 'omega', 'ma'])

# print(data['a'])

# ivp = InvariantPlane(4, data['mass'], data['a'],
#                     data['e'], data['inc'], data['Omega'], data['omega'])

# print(np.rad2deg(ivp.I_Omega_forced_da(10000, 1e-3, np.deg2rad(20), 0)))

# alist = np.linspace(30, 60, 300)
# qlist = []
# plist = []
# for a in alist:
#     q0, p0 = ivp.q_p_forced_da(0, a, 1e-3, np.deg2rad(20), 0)
#     qlist.append(q0)
#     plist.append(p0)
#     print(q0, p0)
# qlist = np.array(qlist)
# plist = np.array(plist)

# fig, ax1 = plt.subplots()
# ax1.scatter(alist, np.rad2deg(qlist))
# ax1.scatter(alist, np.rad2deg(plist))
# ax1.set_ylim(-6, 6)
# ax1.set_xlim(30, 60)
# fig.savefig("I_Omega_forced_da_inc20e0.png", dpi=200)