import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
import os
import pandas as pd
import matplotlib
import orbitplotkit as opk

font = {'weight': 'bold',
        'size': 14}
matplotlib.rc('font', **font)
light_blue = '#5ed2f1'
dark_blue = '#0084c6'
yellow = '#ffb545'
red = '#ff2521'
aN = 30.11

# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT1-Migration-2stage-2EM-400au-2/"
# folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUN-Giant-Synthetic-Ref/"
output_folder = folder + "pics/histogram/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


def plotHistogram(Time):

    # idx = 226
    # time = [0]
    par_size = 0.25
    alpha = 1
    aleft = 20
    a0 = 200
    a1 = 1200
    q0 = 21
    q1 = 105
    # q1 = 600
    qcut = 38
    incmax = 75
    for t in time:
        tri_txt = "reoutput/hist/pl_5.txt"
        df_tri = opk.hist2df(folder+tri_txt)



        pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
        df_pa = opk.snapshot2df(folder+pa_txt)

        fig, ax1 = plt.subplots(figsize=(7, 5))
        print(df_pa['q'].max())


        df_highq = (df_pa['a'] > 50) & opk.isDetached(
            df_pa['a'], df_pa['q'], aN=aN, qcut=qcut)
        df_det = df_pa[df_highq]

        nbins = 11*2
        inner, outer = 100, 1200
        bin_width = (outer - inner)/nbins
        ax1.hist(df_det['a'], nbins, density=True, histtype='stepfilled',facecolor=opk.LIGHT_BLUE,
                           cumulative=False, range=(inner,outer),label='Detached distribution')
        ax1.hist(df_tri['a'], nbins, density=True, histtype='step',edgecolor=opk.DARK_RED,
                           cumulative=False, range=(inner,outer),label='Rogue history')

        # ax1.set_xlim(inner,outer)
        ax1.set_ylabel("Percentage")
        ax1.set_xlabel("a (au)")
        ax1.legend(loc='upper right',fontsize=12)

        yvals = ax1.get_yticks()
        ax1.set_yticklabels(["{0:.0%}".format(y*bin_width) for y in yvals])

        fig.tight_layout()
        plt.savefig(output_folder +
                    "histogram.jpg", dpi=500)
        print("Saved! frame: histogram")
        plt.close()


idx = 621
time_step = 100000000
t0 = int(idx*time_step)
# idx = 31500000000/time_step
# plotSnapshotAQI(time, idx, plotReal=False, plotCold=400000)
time = [t0]
# idx = time[0]/time_step
plotHistogram(time)