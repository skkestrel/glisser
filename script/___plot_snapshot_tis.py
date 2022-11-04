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
output_folder = folder + "pics/tis/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


def plot(time, idx, plotReal=False, plotCold=400000):

    # idx = 226
    # time = [0]
    aleft = 22
    a0 = 200
    a1 = 1200
    T0 = 1.0
    T1 = 5.0
    qcut = 38
    incmax = 55
    for t in time:
        pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
        df_pl = opk.snapshot2df(folder+pl_txt)
        # df_pl['Omega'] = np.deg2rad(df_pl['Omega'])
        pl_num = len(df_pl.index)
        df_pl['color'] = [opk.PINK, opk.ORANGE,
                          opk.DARK_GREEN, opk.RED, opk.ORANGE, opk.LIGHT_BLUE][:pl_num]
        a_n = df_pl['a'][4]

        pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
        df_pa = opk.snapshot2df(folder+pa_txt, apl = a_n)


        fig = plt.figure(figsize=(8.5, 6.5))
        fig.suptitle("Time: {time:6.3f} Myr".format(
            time=t/365.25/1e6), fontsize=15, y=0.95)

        gs = GridSpec(1, 2, width_ratios=[1, 1.6])
        ax3 = fig.add_subplot(gs[0])
        ax1 = fig.add_subplot(gs[1])
        # ax4 = fig.add_subplot(gs[2], sharex=ax3)
        # ax2 = fig.add_subplot(gs[3], sharex=ax1)
        # fig, [ax3, ax4, ax1, ax2] = plt.subplots(2, 2, figsize=(8.5, 6), gridspec_kw={
        #     'height_ratios': [3, 1], 'width_ratios': [2, 1]}, sharex=True)

        # ax3.scatter(df_pl['a'], df_pl['Omega'], alpha=1,
        #             s=150, color=df_pl['color'], marker='x', lw=3, zorder=2)

        par_size = 0.5
        alpha = 1
        df_pa['size'] = par_size
        df_pa['color'] = opk.DARK_GRAY
        df_pa.loc[df_pa.index > plotCold, 'color'] = opk.DARK_RED
        df_highq = (df_pa['a'] > 50) & opk.isDetached(
            df_pa['a'], df_pa['q'], aN=aN, qcut=qcut)

        df_pa.loc[df_highq, 'size'] = par_size*10
        df_pa.loc[df_highq, 'color'] = opk.LIGHT_BLUE
        df_pa = df_pa.sort_values(by=['q'])

        # print(df_highq[df_highq.index < 200000])
        # df_highq[df_highq.index > 200000].to_csv("JSUNT1-Migration-1EM-outer-particles.csv")
        # print(df_highq[df_highq.index > 200000])
        acrit = np.linspace(50, a1)
        qcrit = np.clip(np.nan_to_num(opk.qcrit(acrit, aN),
                        nan=0.0), a_min=qcut, a_max=None)


        # klist, jlist, a_rlist = opk.genOuterRes(
        #     a_n, a_n+1, 300, high1=2, high2=10, order_lim=9)
        # for k, j, a_res in zip(klist, jlist, a_rlist):
        #     label = "{0:0d}/{1:0d}".format(j, k)
        #     if a_res < a0:
        #         ax3.text(a_res, T1+2.5, label,
        #                  fontsize=k**(-0.4)*9, rotation=60, ha='center')
        #     else:
        #         ax1.text(a_res, T1+2.5, label,
        #                  fontsize=k**(-0.4)*9, rotation=60, ha='center')
        #     for ax in [ax3, ax1]:
        #         ax.plot([a_res, a_res], [a_n, T1], ls='dashed',
        #                 color=opk.GRAY, zorder=-1, lw=0.5, alpha=0.8)

        for ax in [ax3, ax1]:
            # ax.scatter(df_pl['a'], df_pl['T'], alpha=1, s=100,
            #            color=df_pl['color'], marker='x', lw=3, zorder=20, rasterized=True)
            ax.scatter(df_pa['a'], df_pa['T'], alpha=alpha,
                       s=df_pa['size'], edgecolors='none', facecolors=df_pa['color'], rasterized=True)

        # ax3.scatter(df_pa['a'], df_pa['Omega'], alpha=alpha,
        #             s=par_size, edgecolors='none', facecolors='C0')
        # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
        # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

        x = np.linspace(aleft, a0, 1000)
        Tis = a_n/x + 2*np.sqrt(x/a_n)
        y = np.linspace(aleft, a0*10000, 1000)
        for ax in [ax3, ax1]:
            ax.plot(x, Tis, color=opk.BLACK, linestyle='dashed', lw=1)
            ax.fill_between(x, Tis, y, facecolor=opk.GRAY, zorder=10)

        if plotReal:
            df_real = pd.read_csv('distant_tno_bary.csv')
            df_real['alpha'] = 1
            df_real.loc[(df_real['a_bary'] < 200), 'alpha'] = 0.5
            df_real['size'] = 30
            df_real.loc[(df_real['a_bary'] < 200), 'size'] = 10
            alpha = df_real['a_bary']/a_n
            df_real['T'] = 1/(alpha) + 2*np.sqrt(alpha*(1-df_real['e_bary']**2))*np.cos(df_real['inc_bary'])
            print(df_real['T'])
            df_real = df_real[(df_real['a_bary'] > 50) & opk.isDetached(
                df_real['a_bary'], df_real['q_bary'], aN=aN, qcut=qcut)]
            for ax in [ax3, ax1]:
                ax.scatter(df_real['a_bary'], df_real['T'], alpha=df_real['alpha'], s=df_real['size'],
                           color=opk.DARK_RED, marker='x', lw=1.5, zorder=2)

        # ax1.plot([30, outer], [38, 38], color=red, alpha=0.5, linestyle='dashed',lw=2)
        # ax2.plot([30, outer], [6, 6], color=red, linestyle='dashed',lw=2)
        # ax1.plot([80, outer], [80, 80], color=red, alpha=0.3, linestyle='dashed',lw=2)
        # for pos,label in zip([aN,38], ["q = 30", "q = 38"]):
        #     ax1.text(15, pos-2, label, fontsize = 12, color=red)

        # rect2 = patches.Rectangle((4, 80), 46, 20, edgecolor='none', facecolor='green', alpha=0.2)
        # ax2.add_patch(rect2)

        # ress = [11/5, 7/3, 5/2, 8/3, 11/4, 3/1, 13/4, 7/2, 15/4, 4/1, 9/2, 5/1]
        # labels = ["11/5", "7/3", "5/2", "8/3", "11/4", "3/1", "13/4","7/2", "15/4", "4/1","9/2", "5/1"]
        # for res, label in zip(ress, labels):
        #     loc = aN*(res)**(2/3)
        #     ax1.plot([loc, loc], [0, loc], color='grey', alpha=0.4 ,linestyle='dashed',lw=1, zorder=-1)
        #     ax2.plot([loc, loc], [0, 120], color='grey', alpha=0.4 ,linestyle='dashed',lw=1, zorder=-1)
        #     ax1.text(loc-0.6, 66.5, label, fontsize = 12, rotation=60)
        ax3.set_axisbelow(False)
        ax1.set_xlabel("a (au)")
        ax3.set_ylabel("Tisserand w/ Neptune (T_n)")
        # ax3.set_ylabel("Omega (deg)")

        ax1.set_xscale("log")
        # ax1.set_yscale("log")
        # ax3.set_yscale("log")
        ax1.set_xlim(a0, a1)

        ax3.set_ylim(T0, T1)
        ax1.set_ylim(T0, T1)
        ax3.set_xlim(aleft, a0)

        # yticks = np.arange(20, round(q1/10-1)*10+1, 10)
        # yticks[0] = 25
        # for ax in [ax1, ax3]:
        #     ax.set_yticks(yticks)
        #     ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        xticks = np.arange(200, 700+1, 100)
        xticks[-1] = 1000
        ax1.set_xticks(xticks)
        ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        # ax1.axes.get_yaxis().set_visible(False)
        ax1.yaxis.tick_right()
        ax1.axes.get_yaxis().set_ticklabels([])

        # ax1.tick_params(labelleft=False, length=0)
        # ax2.tick_params(labelleft=False, length=0)
        # ax2.tick_params(axis="y", labelleft=False)

        # for ax in [ax1, ax2]:
        #     ax.spines['left'].set_visible(False)
        for ax in [ax3]:
            ax.spines['right'].set_visible(False)
        # ax2.set_yscale("log")
        # ax3.set_xlim(a0, a1)
        # ax3.set_ylim(-180, 180)

        fig.tight_layout()
        ax1.xaxis.set_label_coords(0.1, -0.08)
        fig.subplots_adjust(hspace=0, wspace=0)
        plt.savefig(output_folder +
                    "frame_{idx:04d}.jpg".format(idx=idx), dpi=500)
        print("Saved! frame: {idx:04d}".format(idx=idx))
        plt.close()
        idx += 1


idx = 621
time_step = 100000000
t0 = int(idx*time_step)
time = np.arange(t0, 2*36000000001, time_step)
time = [t0]
plot(time, idx, plotReal=True, plotCold=400000)
