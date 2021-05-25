# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.ticker as ticker


# plt.style.use('aiur')

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

COLORS = [BLUE, RED, GREEN, YELLOW]

# ------------------------------
# Semi-major Axes of Nine Planets (PLUTO IS A PLANET!)
# Ratio to Earth Values (or in AU)
A_MERCURY = 0.387
A_VENUS = 0.723
A_EARTH = 1
A_MARS = 1.524
A_JUPITER = 5.20
A_SATURN = 9.58
A_URANUS = 19.20
A_NEPTUNE = 30.05
A_PLUTO = 39.48
A_PLANETS = [A_MERCURY, A_VENUS, A_EARTH, A_MARS,
             A_JUPITER, A_SATURN, A_URANUS, A_NEPTUNE]

# Mass of Nine Planets
# Ratio to Earth Values
M_SUN = 333060.402
M_MERCURY = 0.0553
M_VENUS = 0.815
M_EARTH = 1
M_MARS = 0.107
M_JUPITER = 317.8
M_SATURN = 95.2
M_URANUS = 14.5
M_NEPTUNE = 17.1
M_PLUTO = 0.0025
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

# Diameters of Nine Planets
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

# Diameters of Hill Sphere of Nine Planets
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

# Handy constants
G = 6.67408e-11
YEAR_IN_SEC = 365.25 * 24 * 3600
YEAR_IN_DAY = 365.25
CENT_IN_SEC = 100*YEAR_IN_SEC
CENT_IN_DAY = 100*YEAR_IN_DAY
AU = 149597870700
ARCSEC_CIRCLE = 1296000
RAD_TO_ARCSEC = 206264.806

# ------------------------------


def meanMotion(M, a, G=1):
    return np.sqrt(G*M/a**3)

# ------------------------------


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def fmt(x, pos):
    x *= 100
    return "{0:0.4f}%".format(x)


def wrapTo360(phi):
    return phi % 360


def resonanceLocation(planet, p, q):
    try:
        a = A_PLANETS[planet]
        return (q/p) ** (2/3) * a
    except IndexError:
        print("ERROR!! Planet index number must be in the range of (0,8)")


def resonantAngle(ome1, Ome1, M1, ome2, Ome2, M2, kj, k):
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
