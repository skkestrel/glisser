import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors

def dense_scatter(ax, x, y, data, ny=400, nx=400, logBar=False, label=None, order=None, defaultValue=np.nan, defaultColor=None, upperLimitColor=None, use_cbar=True, log_y=False, ymin=None, vmax=None):
    xmax = x.max()
    xmin = x.min()
    ymax = y.max()

    if not ymin:
        ymin = y.min()

    if log_y:
        yedges, xedges = np.exp(np.linspace(np.log(ymin), np.log(ymax), ny+1)), np.linspace(xmin, xmax, nx+1)
    else:
        yedges, xedges = np.linspace(ymin, ymax, ny+1), np.linspace(xmin, xmax, nx+1)

    H = np.empty((nx, ny))
    H[:] = defaultValue
    X, Y = np.meshgrid(xedges, yedges)

    if log_y:
        dx = (xmax - xmin) / nx
        dy = (np.log(ymax) - np.log(ymin)) / ny
    else:
        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny

    def load(x, y):
        if log_y:
            xc = int(math.floor((xi - xmin) / dx))
            yc = int(math.floor((np.log(yi) - np.log(ymin)) / dy))
            if yc == ny: yc -= 1
            if xc == nx: xc -= 1
        else:
            xc = int(math.floor((xi - xmin) / dx))
            yc = int(math.floor((yi - ymin) / dy))
            if yc == ny: yc -= 1
            if xc == nx: xc -= 1
        return xc, yc

    if order is not None:
        for i in order:
            xi = x[i]
            yi = y[i]
            dat = data[i]
            tup = load(xi, yi)
            if tup[1] < 0:
                continue
            H[tup[0], tup[1]] = dat
    else:
        for xi, yi, dat in zip(x, y, data):
            tup = load(xi, yi)
            if tup[1] < 0:
                continue
            H[tup[0], tup[1]] = dat


    H = H.T[::-1, :]
    ax.grid(False)

    if not vmax:
        vmax = data.max()

    import copy
    cm = copy.copy(matplotlib.cm.get_cmap("viridis"))

    if defaultColor:
        cm.set_bad(color=defaultColor)

    if not use_cbar:
        cmap = None

    if upperLimitColor is not None:
        vmax = np.nextafter(vmax, np.NINF)
        # vmax *= 0.99
        cm.set_over(color=upperLimitColor)

    if logBar:
        norm = matplotlib.colors.LogNorm(vmax=vmax)
    else:
        norm = matplotlib.colors.Normalize(vmax=vmax, clip=upperLimitColor is None)

    if log_y:
        X, Y = np.meshgrid(xedges, yedges)
        H = np.ma.masked_array(H, H < 0)
        a = ax.pcolormesh(X, Y, H[::-1], zorder=0, cmap=cm)
        ax.set_yscale("log")
    else:
        if np.isnan(defaultValue):
            H = np.ma.masked_array(H, np.isnan(H))
        a = ax.imshow(H, aspect='auto', zorder=0, extent=(xmin, xmax, ymin, ymax), norm=norm, cmap=cm)

    if use_cbar:
        cbar = plt.colorbar(a, ax=ax)
        cbar.set_label(label)

    return a

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
