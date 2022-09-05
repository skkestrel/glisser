# -*- coding: UTF-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import matplotlib.font_manager as font_manager
from math import gcd as bltin_gcd
from astropy.io import fits
from astropy.table import Table

plt.style.use('aiur')
# the location of the font file
font_path = '/home/yhuang/GLISSER/roboto/PingFang.ttc'
font_manager.fontManager.addfont(font_path)
roboto = font_manager.FontProperties(fname=font_path)


font = {'family': roboto.get_name(),
        'size': 16}
mpl.rc('font', **font)
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = roboto.get_name()
# mpl.rcParams['mathtext.rm'] = roboto.get_name()

# Pre-defined string
STRING_SMA = 'Semimajor axis (au)'
STRING_INC = 'I (deg)'
STRING_ECC = 'e'
STRING_TIME = 'Time (year)'

# Constants for plotting
DEFAULT_CMAP = "RdYlBu"
SECONDARY_CMAP = "coolwarm_r"
DEFAULT_FONT = "Roboto"
CMAP_START = 0.1
CMAP_END = 1.0
CMAP_BIT = 100
DEFAULT_MARKER_SIZE = 1.5
DEFAULT_LINE_WIDTH = 1.5
DEFAULT_ALPHA = 1

# Custom color palettes
YELLOW = "#FFE026"  # Khala Yellow
BLUE = "#3390B5"  # Khala Blue
LIGHT_BLUE = "#51B9DE"  # Khala Light Blue
GREEN = "#77CC68"  # Aiur Green
RED = "#E72C2C"  # Aiur Red
DARK_RED = "#db5d4d"
LIGHT_YELLOW = "#fdecbc"
DARK_GREEN = "#439e64"
ORANGE = "#ffb545"
PINK = "#f10075"
GRAY = "#AEAEAE"
DARK_GRAY = "#585858"
BLACK = "#3B3B3B"
GREY = GRAY

COLORS = [BLUE, RED, GREEN, YELLOW]

# ------------------------------
# Semi-major Axes of 8 Planets (PLUTO IS NOT A PLANET!)
# Ratio to Earth Values (or in AU)
A_MERCURY = 0.387
A_VENUS = 0.723
A_EARTH = 1
A_MARS = 1.524
A_JUPITER = 5.20
A_SATURN = 9.54
A_URANUS = 19.19
A_NEPTUNE = 30.07
A_PLUTO = 39.48
A_PLANETS = [A_MERCURY, A_VENUS, A_EARTH, A_MARS,
             A_JUPITER, A_SATURN, A_URANUS, A_NEPTUNE]

# Mass of 8 Planets
# Ratio to Earth Values
M_SUN = 1.0
M_MERCURY = 1.6601141530543488e-07
M_VENUS = 2.4478382877847715e-06
M_EARTH = 3.040432648022642e-06
M_MARS = 3.2271560375549977e-07
M_JUPITER = 0.0009547919152112404
M_SATURN = 0.0002858856727222417
M_URANUS = 4.36624373583127e-05
M_NEPTUNE = 5.151383772628674e-05
M_PLANETS = [M_MERCURY, M_VENUS, M_EARTH, M_MARS,
             M_JUPITER, M_SATURN, M_URANUS, M_NEPTUNE]

# Mass of the solar system bodies
# SI units
m_sun = 1.9890E+30
m_mercury = 3.3011E+23
m_venus = 4.8675E+24
m_earth = 5.9724E+24
m_mars = 6.4171E+23
m_jupiter = 1.8982E+27
m_saturn = 5.6834E+26
m_uranus = 8.6813E+25
m_neptune = 1.0241E+26
m_planets = [m_mercury, m_venus, m_earth, m_mars,
             m_jupiter, m_saturn, m_uranus, m_neptune]

# Inclinations of 8 Planets
# units: degrees
I_MERCURY = 7.00
I_VENUS = 3.39
I_EARTH = 0.0
I_MARS = 1.85
I_JUPITER = 1.31
I_SATURN = 2.48
I_URANUS = 0.77
I_NEPTUNE = 1.77
I_PLUTO = 17.14
I_PLANETS = [I_MERCURY, I_VENUS, I_EARTH, I_MARS,
             I_JUPITER, I_SATURN, I_URANUS, I_NEPTUNE]


# Diameters of 8 Planets
# Ratio to Earth Values
D_MERCURY = 0.383
D_VENUS = 0.949
D_EARTH = 1
D_MARS = 0.532
D_JUPITER = 11.21
D_SATURN = 9.45
D_URANUS = 4.01
D_NEPTUNE = 3.88
D_PLUTO = 0.186
D_PLANETS = [D_MERCURY, D_VENUS, D_EARTH, D_MARS,
             D_JUPITER, D_SATURN, D_URANUS, D_NEPTUNE]

# Diameters of Hill Sphere of 8 Planets
# in AU
H_MERCURY = 0.0011718
H_VENUS = 0.006713
H_EARTH = 0.009837
H_MARS = 0.00657
H_JUPITER = 0.33805
H_SATURN = 0.412
H_URANUS = 0.44645
H_NEPTUNE = 0.7691
H_PLUTO = 0.03854
H_PLANETS = [H_MERCURY, H_VENUS, H_EARTH, H_MARS,
             H_JUPITER, H_SATURN, H_URANUS, H_NEPTUNE]


# Eigen frenquencies ans the associated phases of the Solay System
# frenquencies in arcseconds per year; phases in degrees
# (Laskar 1990)

# g_freqs = [5.59, 7.455, 17.30, 17.85,
#            4.24882, 28.2203, 3.08952, 0.66698]
# f_freqs = [-5.59, -7.00, -18.88, -17.80,
#            0.0, -26.33020, -3.00563, -0.69195]
g_freqs = [5.46, 7.34, 17.32, 18.00,
           4.30, 27.77, 2.72, 0.63332]
f_freqs = [-5.20, -6.57, -18.74, -17.63,
           0.0, -25.73, -2.90, -0.67752]
beta_phases = [92.182, 196.881, 335.224,
               317.948, 29.550, 125.120, 131.944, 69.021]
gamma_phases = [19.433, 318.057, 255.031,
                296.514, 107.102, 127.367, 315.063, 202.293]

g1 = 5.59
g2 = 7.453
g3 = 17.368
g4 = 17.916
g5 = 4.257482
g6 = 28.2449
g7 = 3.087946
g8 = 0.673019

s1 = -5.61
s2 = -7.06
s3 = -18.848
s4 = -17.751
s5 = 0.0
s6 = -26.347841
s7 = -2.9925258
s8 = -0.691740


Iij = [12449, 1178, 849, 180, -2, -3, 2, 0,
       -3545, 1004, 810, 180, -1, -2, 1, 0,
       409, -2680, 2448, -3589, 0, 0, 0, 0,
       116,  -685,  453, 5025, 0, -2, 0, 0,
       2757, 2757, 2757, 2757, 2757, 2757, 2757, 2757,
       28, 12, 281, 965, -631, 1572, -69, -8,
       -332, -192, -173, -126, -96, -78, 1760, -207,
       -145, -132, -130, -123, -117, -113, 110, 1175]

Iij = np.array(Iij).reshape(8, 8)*1e-5
Iji = Iij.T

# Handy constants
G = 6.67408e-11
YEAR_IN_SEC = 365.25 * 24 * 3600
YEAR_IN_DAY = 365.25
CENT_IN_SEC = 100*YEAR_IN_SEC
CENT_IN_DAY = 100*YEAR_IN_DAY
AU = 149597870700
ARCSEC_CIRCLE = 1296000
RAD_TO_ARCSEC = 206264.806
ARCSEC_TO_RAD = 1/RAD_TO_ARCSEC

# ------------------------------


def qcrit(a, aN):
    return aN*np.sqrt(np.log(115.2*m_neptune/m_sun*(a/aN)**(2.5)))


def isDetached(a, q, aN, qcut=38):
    qc = np.nan_to_num(qcrit(a, aN), nan=0.0)
    qc = np.clip(qc, a_min=qcut, a_max=None)
    return q > qc


def snapshot2df(txtfile, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' '):
    df = pd.read_csv(txtfile, names=names, delimiter=delimiter).set_index('id')
    df['q'] = df['a'] * (1 - df['e'])
    df['ascnode'] = (df['a'] * (1 - df['e']**2)) / \
        (1+df['e']*np.cos(-df['omega']))
    df['inc'] = np.rad2deg(df['inc'])
    df['Omega'] = np.rad2deg(df['Omega'])
    df['omega'] = np.rad2deg(df['omega'])
    df['f'] = np.rad2deg(df['f'])
    return df


def fits2csv(fitsfile, output):
    with fits.open(fitsfile, memmap=True) as hdu:

        print(hdu)
        # read into an astropy Table object
        table = Table(hdu[1].data)

        # write to a CSV file
        table.write(output, delimiter=',', format='ascii', overwrite=True)


def mpc2Compact(code):
    code = code.strip('C/')
    code = code.replace(" ", "")
    year = int(code[:2])
    if year == 18:
        first = 'I'
    elif year == 19:
        first = 'J'
    elif year == 20:
        first = 'K'

    year2 = code[2:4]
    first_letter = code[4]
    second_letter = code[5]
    if len(code) == 6:
        middle = '00'
    elif len(code) == 7:
        middle = '0' + code[6]
    elif len(code) == 8:
        middle = code[6:]
    elif len(code) == 9:
        num1 = int(code[6:8])
        if num1 < 36:
            letter = chr(num1 + 55)
        else:
            letter = chr(num1 + 61)
        middle = letter + code[8]

    return first+year2+first_letter+middle+second_letter


def Compact2MPC(code):
    year = ord(code[0]) - 55
    if year < 18 or year > 20 or len(code) != 7:
        return '-'
    year = int(code[1:3]) + year * 100
    first_letter, second_letter, last_digit, double_digit = code[3], code[6], int(
        code[5]), code[4]
    if double_digit.isdigit():
        number = 10*int(double_digit) + last_digit
        if number == 0:
            number = ''
    elif double_digit.isupper():
        number = 10*(ord(double_digit) - 55) + last_digit
    else:
        number = 10*(ord(double_digit) - 61) + last_digit
    return (str(year) + ' ' + first_letter + second_letter + str(number))

# ------------------------------


def aProbFuncNorm(a, arange, aco):
    x1, x2, co, xx = arange[0], arange[1], aco, a
    scale = -1/(1+co) * x1**(1+co) + 1/(1+co) * x2**(1+co)
    # scale = x1**co
    # print(scale)
    return xx**(co)/scale


def coprime(a, b):
    return bltin_gcd(a, b) == 1


def meanMotion(M, a, G=1):
    return np.sqrt(G*M/a**3)


def F2M(f, e):
    tanE2 = np.tan(f/2)*np.sqrt((1-e)/(1+e))
    E2 = np.arctan(tanE2)
    E = E2 * 2
    return E2M(E, e)


def E2M(E, e):
    return E - e*np.sin(E)


def elements2R(a, e, f):
    return a*(1-e**2)/(1+e*np.cos(f))


def M2F(M, e):
    # input: radians
    # output: radians
    acc = 1e-14
    if (M == 0.0):
        return 0.0
    if (M < np.pi and M > 0):
        E = M + e / 2
    else:
        E = M - e / 2
    incr = 1
    count = 0
    while (np.abs(incr) > acc):
        incr = keplerIteration(E, e, M)
        E += incr
    result = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))
    return result


def keplerIteration(E, e, M):
    return -(E - e * np.sin(E) - M) / (1 - e * np.cos(E))


def vinf(a, e, inc, apl):
    a = a/apl
    vinf2 = 3 - 1/a - 2*np.sqrt(a*(1-e**2))*np.cos(inc)
    return np.sqrt(vinf2)


def vpa(a, apl):
    a = a/apl
    return np.sqrt(2 - 1/a)


def alpha(vinf, vpa):
    return np.arccos((vpa**2-vinf**2-1)/(2*vinf))

# ------------------------------


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def percent_fmt(x, pos):
    x *= 100
    return "{0:0.0f}%".format(x)


def wrapTo360(phi):
    return phi % 360


def resonanceLocation(planet, p, q):
    try:
        a = A_PLANETS[planet]
        return (q/p) ** (2/3) * a
    except IndexError:
        print("ERROR!! Planet index number must be in the range of (0,8)")


def resonanceLocationBYA(a, p, q):
    return (q/p) ** (2/3) * a


def genOuterRes(a_pl, a1, a2, high1=4, high2=30, order_lim=15):
    klist, jlist = [], []
    a_rlist = []
    for h1 in range(1, high1+1):
        for h2 in range(1, high2+1):
            if coprime(h1, h2) and abs(h1 - h2) <= order_lim:
                a_res = resonanceLocationBYA(a_pl, h1, h2)
                if a_res < a2 and a_res > a1:
                    klist.append(h1)
                    jlist.append(h2)
                    a_rlist.append(a_res)
    return klist, jlist, a_rlist


def resonantAngle(Ome1, ome1, M1, Ome2, ome2, M2, kj, k):
    varpi1 = ome1 + Ome1
    varpi2 = ome2 + Ome2
    lambda1 = varpi1 + M1
    lambda2 = varpi2 + M2
    order = np.abs(kj-k)
    if order == 0:
        phi = lambda1 - lambda2
    elif k > 0:
        phi = kj * lambda2 - k * lambda1 - (kj - k) * varpi1
        phi /= kj-k
    elif k < 0:
        k = np.abs(k)
        phi = k * lambda1 - kj * lambda2 - (kj + k) * varpi1
        phi /= kj+k
    return phi


def resonantAngleOuter(O_pl, o_pl, M_pl, O_pa, o_pa, M_pa, j, k):
    # for a j:k outer resonance (j > k)
    varpi_pl = O_pl + o_pl
    varpi_pa = O_pa + o_pa
    lambda_pl = varpi_pl + M_pl
    lambda_pa = varpi_pa + M_pa
    phi_jk = j * lambda_pa - k * lambda_pl - \
        (j-k) * varpi_pa  # Gladman et al. (2012)
    return phi_jk

# ------------------------------


def scatterPlot(ax, x, y, xlable, ylabel, title, xlim, ylim):
    sp = ax.scatter(x, y)
    setBasicLabels(ax, xlable, ylabel, title, xlim, ylim)
    return sp


def scatterPlotWithSize(ax, x, y, sz, xlable, ylabel, title, xlim, ylim):
    sp = ax.scatter(x, y, s=sz, edgecolors='black', linewidth=0.5)
    setBasicLabels(ax, xlable, ylabel, title, xlim, ylim)
    return sp


def scatterPlotWithSizeAndColor(ax, x, y, sz, c, xlable, ylabel, title, xlim, ylim):
    sp = ax.scatter(x, y, s=sz, c=c, cmap=SECONDARY_CMAP, norm=colors.LogNorm(
        vmin=c.min(), vmax=c.max()), edgecolors='black', linewidth=0.5)
    setBasicLabels(ax, xlable, ylabel, title, xlim, ylim)
    return sp


def scatterColorPlot(ax, x, y, c, xlable, ylabel, title, xlim, ylim):
    scp = ax.scatter(x, y, c=c, cmap=DEFAULT_CMAP,
                     norm=mpl.colors.LogNorm(vmin=1e6, vmax=2.5e9))
    setBasicLabels(ax, xlable, ylabel, title, xlim, ylim)
    return scp


def setBasicLabels(ax, xlable, ylabel, title, xlim, ylim):
    ax.set_xlabel(xlable)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])


def resonancePlot(ax, a_res, ylim, color, alpha):
    rp = ax.plot((a_res, a_res), ylim, color=color,
                 linestyle='dashed', alpha=alpha, zorder=-4)
    return rp


def stdTwoScatterPlot(axes, alltp, alim, elim, Ilim):
    a, e, I = alltp[:, 0], alltp[:, 1], alltp[:, 2]

    scatterPlot(axes[0], a, e, r'', r'e', '', alim, elim)
    scatterPlot(axes[1], a, I, r'Semi-major axis (AU)',
                r'I (deg)', '', alim, Ilim)


def stdFourScatterPlot(axes, alltp, alim, elim, Ilim):
    a, e, I, ome, Ome = alltp[:, 0], alltp[:,
                                           1], alltp[:, 2], alltp[:, 3], alltp[:, 4]

    scatterPlot(axes[0, 0], a, e, 'a', 'e', 'a-e', alim, elim)
    scatterPlot(axes[0, 1], a, I, 'a', 'I', 'a-I', alim, Ilim)
    scatterPlot(axes[1, 0], e, ome, 'e', 'omega',
                'e-omega', elim, (0, 360))
    scatterPlot(axes[1, 1], I, Ome, 'I', 'Omega',
                'I-Omega', Ilim, (0, 360))


def drawABox(ax, xx, yy, color='pink', alpha=0.5):
    left, bottom, width, height = (xx[0], yy[0], xx[1]-xx[0], yy[1]-yy[0])
    rect = plt.Rectangle((left, bottom), width, height, linewidth=DEFAULT_LINE_WIDTH, facecolor=color,
                         fill=True, alpha=alpha, zorder=-1)
    ax.add_patch(rect)
    return rect


def drawABoxOutlines(ax, xx, yy, color='pink', alpha=0.5):
    left, bottom, width, height = (xx[0], yy[0], xx[1]-xx[0], yy[1]-yy[0])
    rect = plt.Rectangle((left, bottom), width, height, linewidth=DEFAULT_LINE_WIDTH, facecolor='none', edgecolor=color,
                         fill=False, alpha=alpha, zorder=-1)
    ax.add_patch(rect)
    return rect


def NEOSSatModelAEPlot(ax, grid, alim, elim):
    a_step, e_step, I_step = 0.05, 0.02, 2

    a_count = int((alim[1] - alim[0]) / a_step) + 1
    a_start = int((alim[0]-0.025)/a_step) + 1
    a_range = np.linspace(alim[0], alim[1], a_count)
    a_left = alim[0]-a_step/2
    a_right = alim[1]+a_step/2

    e_count = int((elim[1] - elim[0]) / e_step) + 1
    e_start = int((elim[0]-0.01)/e_step)
    e_range = np.linspace(elim[0], elim[1], e_count)
    e_left = elim[0]-e_step/2
    e_right = elim[1]+e_step/2

    new_cmap = truncate_colormap(plt.get_cmap(
        'gist_yarg'), CMAP_START, CMAP_END,  CMAP_BIT)
    heatmap = ax.imshow(grid, extent=[
                        a_left, a_right, e_left, e_right], cmap=new_cmap, zorder=-10)
    cbar = ax.figure.colorbar(heatmap, ax=ax, pad=0.01,
                              format=ticker.FuncFormatter(fmt))
    # ax.set_xlabel("Semi-major axis")
    ax.set_ylabel("e")
    return heatmap


def NEOSSatModelAIPlot(ax, grid, alim, Ilim):
    a_step, e_step, I_step = 0.05, 0.02, 2

    a_count = int((alim[1] - alim[0]) / a_step) + 1
    a_start = int((alim[0]-0.025)/a_step) + 1
    a_range = np.linspace(alim[0], alim[1], a_count)
    a_left = alim[0]-a_step/2
    a_right = alim[1]+a_step/2

    I_count = int((Ilim[1] - Ilim[0]) / I_step) + 1
    I_start = int((Ilim[0]-0.01)/I_step)
    I_range = np.linspace(Ilim[0], Ilim[1], I_count)
    I_left = Ilim[0]-I_step/2
    I_right = Ilim[1]+I_step/2

    new_cmap = truncate_colormap(plt.get_cmap(
        'gist_yarg'), CMAP_START, CMAP_END, CMAP_BIT)
    heatmap = ax.imshow(grid, extent=[
                        a_left, a_right, I_left, I_right], cmap=new_cmap, zorder=-10)
    cbar = ax.figure.colorbar(heatmap, ax=ax, pad=0.01,
                              format=ticker.FuncFormatter(fmt))

    ax.set_xlabel(STRING_SMA)
    ax.set_ylabel(STRING_INC)
    return heatmap


def crossingLinePlot(ax, planetcode, color, alpha=1):
    e = np.linspace(0, 1, 1000)
    ap = A_PLANETS[planetcode]
    a1 = ap/(1+e)
    a2 = ap/(1-e)
    rp = ax.plot(a1, e, color=color, alpha=alpha)
    rp = ax.plot(a2, e, color=color, alpha=alpha)
    return rp


def crossingBandPlotR(ax, arange, color, alpha=0.35):
    # a1 = np.linspace(arange[0], 1.5, 500)
    a1 = np.linspace(1, 1.5, 500)
    e1 = 1 - arange[0]/a1
    e2 = 1 - arange[1]/a1

    rp = ax.fill_between(a1, e1, e2, color=color, alpha=alpha/2, zorder=-5)
    rp = ax.fill_between(a1, e1, 1, color=color, alpha=alpha, zorder=-5)
    rp = ax.plot(a1, e1, color=color, alpha=alpha)
    rp = ax.plot(a1, e2, color=color, alpha=alpha)
    return rp


def crossingBandPlotL(ax, arange, color, alpha=0.35):
    a1 = np.linspace(1, arange[0], 500)
    e1 = arange[0]/a1 - 1
    e2 = arange[1]/a1 - 1

    rp = ax.fill_between(a1, e1, e2, color=color, alpha=alpha/2, zorder=-5)
    rp = ax.fill_between(a1, e1, 1, color=color, alpha=alpha, zorder=-5)
    rp = ax.plot(a1, e1, color=color, alpha=alpha)
    rp = ax.plot(a1, e2, color=color, alpha=alpha)
    return rp
