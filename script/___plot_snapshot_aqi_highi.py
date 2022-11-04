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
# folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT1-Migration-2EM-3/"
folder = "/home/yhuang/GLISSER/"
output_folder = folder
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


def plotSnapshotHighI(time, idx, plotReal=False):

    # idx = 226
    # time = [0]
    par_size = 0.25
    alpha = 1
    aleft = 0.1
    a0 = 100
    a1 = 1200
    q0 = 0.5
    q1 = 50
    qcut = 15
    incmin = 55
    incmax = 185
    for t in time:

        fig = plt.figure(figsize=(8, 5.5))

        gs = GridSpec(2, 2, width_ratios=[1.2, 1], height_ratios=[2, 1])
        ax3 = fig.add_subplot(gs[0])
        ax1 = fig.add_subplot(gs[1])
        ax4 = fig.add_subplot(gs[2], sharex=ax3)
        ax2 = fig.add_subplot(gs[3], sharex=ax1)
 
        x = np.linspace(aleft, max(q1, a0), 1000)
        y = np.linspace(aleft, max(q1, a0)*10000, 1000)
        for ax in [ax3, ax1]:
            ax.plot(x, x, color=opk.BLACK, linestyle='dashed', lw=1)
            ax.fill_between(x, x, y, facecolor=opk.GRAY, zorder=10)

        if plotReal:
            s = 8
            df_real = pd.read_csv('sbdb_query_results_highi_9.22.csv')
            # print(a_n)
            df_real['color'] = opk.GRAY
            df_real['size'] = s
            # df_real.loc[df_real['i'] > 90, 'size'] = s*2
            df_real.loc[df_real['q'] > qcut, 'color'] = opk.DARK_RED
            df_real.loc[df_real['full_name']=="514107 Ka`epaoka`awela (2015 BZ509)", 'color'] = opk.BLUE
            df_real.loc[df_real['q'] > qcut, 'size'] = s*3
            df_real.loc[df_real['full_name']=="514107 Ka`epaoka`awela (2015 BZ509)", 'size'] = s*3
            df_real = df_real.sort_values(by=['size'])
            for ax in [ax3, ax1]:
                ax.scatter(df_real['a'], df_real['q'], alpha=1, s=df_real['size'], edgecolor='none', facecolor=df_real['color'], zorder=2)
                ax.plot([aleft, a1], [qcut, qcut], ls='dashed',
                    color=opk.DARK_RED, zorder=2, alpha=1, lw=1)
            for ax in [ax4, ax2]:
                ax.scatter(df_real['a'], df_real['i'], alpha=1, s=df_real['size'],
                edgecolor='none', facecolor=df_real['color'], zorder=2)


        ax3.set_axisbelow(False)
        ax2.set_xlabel("a (au)")
        ax3.set_ylabel("q (au)")
        ax4.set_ylabel("I (deg)")
        # ax3.set_ylabel("Omega (deg)")

        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax3.set_yscale("log")
        ax2.set_xscale("log")
        ax2.set_ylim(incmin, incmax)
        ax1.set_xlim(a0, a1)

        ax3.set_ylim(q0, q1)
        ax1.set_ylim(q0, q1)
        ax3.set_xlim(aleft, a0)
        # ax2.set_yscale("log")
        ax4.set_ylim(incmin, incmax)

        yticks = [1, 5, 10, 15, 20, 30]
        # yticks[0] = 25
        for ax in [ax1, ax3]:
            ax.set_yticks(yticks)
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        yticks = [60, 90, 120, 150, 180]
        # yticks[0] = 25
        for ax in [ax2, ax4]:
            ax.set_yticks(yticks)
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        xticks = np.arange(200, 601, 100)
        xticks[-2] = 600
        xticks[-1] = 1000
        ax2.set_xticks(xticks)
        ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        # ax1.axes.get_yaxis().set_visible(False)
        ax1.yaxis.tick_right()
        ax1.axes.get_yaxis().set_ticklabels([])
        ax2.yaxis.tick_right()
        ax2.axes.get_yaxis().set_ticklabels([])

        # ax1.tick_params(labelleft=False, length=0)
        # ax2.tick_params(labelleft=False, length=0)
        # ax2.tick_params(axis="y", labelleft=False)

        # for ax in [ax1, ax2]:
        #     ax.spines['left'].set_visible(False)
        for ax in [ax3, ax4]:
            ax.spines['right'].set_visible(False)
        # ax2.set_yscale("log")
        # ax3.set_xlim(a0, a1)
        # ax3.set_ylim(-180, 180)

        fig.align_ylabels([ax3, ax4])
        fig.tight_layout()
        ax2.xaxis.set_label_coords(0, -0.18)
        fig.subplots_adjust(hspace=0, wspace=0)
        plt.savefig(output_folder +
                    "frame_{idx:04d}.jpg".format(idx=idx), dpi=500)
        print("Saved! frame: {idx:04d}".format(idx=idx))
        plt.close()
        idx += 1


idx = 0
time_step = 50000000
t0 = int(idx*time_step)
time = np.arange(t0, 2*37000000001, time_step)

# idx = 31500000000/time_step
# plotSnapshotAQI(time, idx, plotReal=False, plotCold=400000)
time = [t0]
# idx = time[0]/time_step
plotSnapshotHighI(time, idx, plotReal=True)