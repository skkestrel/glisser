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
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT1-Migration-1EM/"
output_folder = folder + "pics/snapshots/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


def plot(time, idx, plotReal=False):

    # idx = 226
    # time = [0]
    aleft = 3
    a0 = 100
    a1 = 1200
    q0 = 3
    q1 = 105
    qcut = 40
    incmax = 75
    for t in time:
        pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
        df_pl = opk.snapshot2df(folder+pl_txt)
        # df_pl['Omega'] = np.deg2rad(df_pl['Omega'])
        pl_num = len(df_pl.index)
        df_pl['color'] = [opk.PINK, opk.ORANGE,
                          opk.DARK_GREEN, opk.RED, opk.LIGHT_BLUE, opk.BLUE][:pl_num]
        pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
        df_pa = opk.snapshot2df(folder+pa_txt)
        # df_pa['Omega'] = np.deg2rad(df_pa['Omega'])

        fig = plt.figure(figsize=(10, 6.5))
        fig.suptitle("Time: {time:6.3f} Myr".format(
            time=t/365.25/1e6), fontsize=15, y=0.95)

        gs = GridSpec(2, 2, width_ratios=[1, 2], height_ratios=[4, 1])
        ax3 = fig.add_subplot(gs[0])
        ax1 = fig.add_subplot(gs[1])
        ax4 = fig.add_subplot(gs[2], sharex=ax3)
        ax2 = fig.add_subplot(gs[3], sharex=ax1)

        fig.subplots_adjust(hspace=0, right=0.8, left=0.08)
        pos1 = ax1.get_position()
        # pos2 = ax4.get_position()
        # print(pos1)
        rect_histy = [pos1.x1, pos1.y0, 0.15, pos1.y1 - pos1.y0]
        ax5 = fig.add_axes(rect_histy, sharey=ax3)
        # fig, [ax3, ax4, ax1, ax2] = plt.subplots(2, 2, figsize=(8.5, 6), gridspec_kw={
        #     'height_ratios': [3, 1], 'width_ratios': [2, 1]}, sharex=True)

        # ax3.scatter(df_pl['a'], df_pl['Omega'], alpha=1,
        #             s=150, color=df_pl['color'], marker='x', lw=3, zorder=2)

        par_size = 0.25
        alpha = 1
        df_pa['size'] = par_size
        df_pa['color'] = opk.DARK_GRAY
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

        a_n = df_pl['a'][4]
        klist, jlist, a_rlist = opk.genOuterRes(
            a_n, a_n+1, 300, high1=2, high2=9, order_lim=8)
        for k, j, a_res in zip(klist, jlist, a_rlist):
            label = "{0:0d}/{1:0d}".format(j, k)
            if a_res < a0:
                ax3.text(a_res, q1+2.5, label,
                        fontsize=k**(-0.4)*9, rotation=60, ha='center')
            else:
                ax1.text(a_res, q1+2.5, label,
                        fontsize=k**(-0.4)*9, rotation=60, ha='center')
            for ax in [ax3, ax1]:
                ax.plot([a_res, a_res], [a_n, q1], ls='dashed',
                        color=opk.GRAY, zorder=-1, lw=0.5, alpha=0.8)

        for ax in [ax3, ax1]:
            ax.scatter(df_pl['a'], df_pl['q'], alpha=1, s=100,
                       color=df_pl['color'], marker='x', lw=3, zorder=20, rasterized=True)
            ax.scatter(df_pa['a'], df_pa['q'], alpha=alpha,
                       s=df_pa['size'], edgecolors='none', facecolors=df_pa['color'], rasterized=True)
            ax.plot(acrit, qcrit, ls='dashed',
                    color=opk.BLUE, zorder=2, alpha=1, lw=1)
            ax.plot([50, 50], [qcut, 50], ls='dashed',
                    color=opk.BLUE, zorder=2, alpha=1, lw=1)
            ax.plot([a_n, a1], [a_n, a_n], ls='dashed',
                    color=opk.DARK_RED, zorder=2, alpha=1, lw=1)

        for ax in [ax4, ax2]:
            ax.scatter(df_pl['a'], df_pl['inc'], alpha=1,
                       s=100, color=df_pl['color'], marker='x', lw=3, zorder=20, rasterized=True)
            ax.scatter(df_pa['a'], df_pa['inc'], alpha=alpha,
                       s=df_pa['size'], edgecolors='none', facecolors=df_pa['color'], zorder=1, rasterized=True)

        ax5.scatter(df_pl['inc'], df_pl['q'], alpha=1,
                       s=100, color=df_pl['color'], marker='x', lw=3, zorder=20, rasterized=True)
        ax5.scatter(df_pa['inc'], df_pa['q'], alpha=alpha,
                       s=df_pa['size'], edgecolors='none', facecolors=df_pa['color'], zorder=1, rasterized=True)

        # ax3.scatter(df_pa['a'], df_pa['Omega'], alpha=alpha,
        #             s=par_size, edgecolors='none', facecolors='C0')
        # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
        # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

        # aN = df_pl['a'].iloc[-2]
        x = np.linspace(aleft, max(q1, a0), 1000)
        y = np.linspace(aleft, max(q1, a0)*10000, 1000)
        for ax in [ax3, ax1]:
            ax.plot(x, x, color=opk.BLACK, linestyle='dashed', lw=1)
            ax.fill_between(x, x, y, facecolor=opk.GRAY, zorder=10)

        if plotReal:
            df_real = pd.read_csv('distant_tno_bary.csv')
            print(a_n)
            # df_real['a_bary']/30.11*a_n
            df_real = df_real[(df_real['a_bary'] > 50) & opk.isDetached(
                df_real['a_bary'], df_real['q_bary'], aN=aN, qcut=qcut)]
            for ax in [ax3, ax1]:
                ax.scatter(df_real['a_bary'], df_real['q_bary'], alpha=1, s=10,
                           color=opk.DARK_RED, marker='x', lw=1, zorder=2)
            for ax in [ax4, ax2]:
                ax.scatter(df_real['a_bary'], np.rad2deg(df_real['inc_bary']), alpha=1, s=10,
                           color=opk.DARK_RED, marker='x', lw=1, zorder=2)
            ax5.scatter(np.rad2deg(df_real['inc_bary']), df_real['q_bary'], alpha=1, s=10,
                           color=opk.DARK_RED, marker='x', lw=1, zorder=2)

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
        ax2.set_xlabel("a (au)")
        ax3.set_ylabel("q (au)")
        ax4.set_ylabel("I (deg)")
        # ax3.set_ylabel("Omega (deg)")

        ax1.set_xscale("log")
        # ax1.set_yscale("log")
        # ax3.set_yscale("log")
        ax2.set_xscale("log")
        ax2.set_ylim(-2, incmax)
        ax1.set_xlim(a0, a1)

        ax3.set_ylim(q0, q1)
        ax1.set_ylim(q0, q1)
        ax3.set_xlim(aleft, a0)
        # ax2.set_yscale("log")
        ax4.set_ylim(-2, incmax)
        ax5.set_xlim(-2, incmax)
        # ax5.set_ylim(4, q1)

        yticks = np.arange(20, round(q1/10-1)*10+1, 10)
        yticks[0] = 25
        for ax in [ax1, ax3]:
            ax.set_yticks(yticks)
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        xticks = np.arange(100, 600+1, 100)
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

        # fig.tight_layout()
        ax2.xaxis.set_label_coords(0.2, -0.35)
        fig.subplots_adjust(hspace=0, wspace=0)
        plt.savefig(output_folder +
                    "frame_qi_{idx:04d}.jpg".format(idx=idx), dpi=500)
        print("Saved! frame: {idx:04d}".format(idx=idx))
        plt.close()
        idx += 1


idx = 1096
time_step = 50000000
t0 = int(idx*time_step)
time = np.arange(t0, 2*36000000001, time_step)
time = [t0]
plot(time, idx, plotReal=True)
