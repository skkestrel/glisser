import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
import struct
import os
import pandas as pd
import matplotlib
import orbitplotkit as opk

font = {'weight': 'bold',
        'size': 18}
matplotlib.rc('font', **font)


# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
# folder1 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNR-Resonant-filter/"
# folder1 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUN-Resonant-filter/"
folder1 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Giant-Synthetic/"
output_folder = folder1 + "/pics/snapshots_RESO"
output_folder = "/home/yhuang/GLISSER/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
aN_bary = 30.068010875914965


def plot(time, plotREAL):
    # df_temp = pd.read_csv(
    #     folder1+"classes.txt", names=['id', 'RESO', 'k1', 'k2'], delimiter=' ').set_index('id')

    # time = [0, 21520000000]
    idx = 1
    # time = [0]

    par_size = 3
    qcut = 38
    q_low_bound = 32.5
    q_up_bound = 69
    i0, i1 = 0, 62
    per1 = 0.12
    outer = 101
    inner = 49.5

    nbin = 200

    for t in time:
        fig, [ax3, ax1, ax2] = plt.subplots(3, figsize=(
            13, 9), gridspec_kw={'height_ratios': [1, 4, 1]}, sharex=True)

        pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
        df_pl = pd.read_csv(
            folder1+pl_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')
        df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
        df_pl['inc'] = np.rad2deg(df_pl['inc'])
        df_pl['Omega'] = np.rad2deg(df_pl['Omega'])
        df_pl['omega'] = np.rad2deg(df_pl['omega'])
        df_pl['M'] = np.rad2deg(df_pl['M'])

        pa_txt = "reoutput/snapshots_RESO/temp_{0:d}.csv".format(t)
        df_pa = pd.read_csv(folder1+pa_txt, delimiter=',').set_index('id')
        df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])
        # df_pa['q0'] = df_pa0['q']
        # df_pa['inc0'] = df_pa0['inc']
        df_pa['inc'] = np.rad2deg(df_pa['inc'])
        df_pa['Omega'] = np.rad2deg(df_pa['Omega'])
        df_pa['omega'] = np.rad2deg(df_pa['omega'])
        df_pa['M'] = np.rad2deg(df_pa['M'])
        # print(df_pa)

        df_pa['alpha'] = 0.8
        # df_pa.loc[df_pa['RESO'] == 1, 'alpha'] = 0.5
        df_pa.loc[df_pa['q'] < qcut, 'alpha'] = 0.25
        # df_pa.loc[df_pa['q'] > 40, 'alpha'] = 1

        df_pa['parsize'] = 1.5
        df_pa.loc[df_pa['q'] > qcut, 'parsize'] = 6
        # df_pa.loc[df_pa['q'] > 40, 'parsize'] = 4

        df_pa['color'] = opk.GRAY

        df_pa.loc[df_pa['q'] > qcut, 'color'] = opk.BLUE
        df_pa.loc[df_pa['RESO'] == -1, 'color'] = opk.DARK_RED
        df_pa.loc[df_pa['RESO'] > 1, 'color'] = opk.ORANGE
        df_pa.loc[df_pa['q'] < qcut, 'color'] = opk.GRAY
        df_pa.loc[df_pa['RESO'] == -2, 'color'] = opk.GRAY

        df_highq = df_pa[df_pa['q'] > qcut].copy()
        df_arange = df_pa[df_pa['a'].between(50, 100)].copy()
        # df_highq = df_highq[df_highq['a'].between(50, 100)]
        # print(df_highq)

        df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
        ax1.scatter(df_pl['a'], df_pl['q'], alpha=0.8, s=100,
                    color=opk.DARK_RED, marker='x', lw=2, zorder=2)
        ax2.scatter(df_pl['a'], np.rad2deg(df_pl['inc']), alpha=0.8,
                    s=100, color=opk.DARK_RED, marker='x', lw=2, zorder=2)

        ax1.scatter(df_pa['a'], df_pa['q'], alpha=df_pa['alpha'],
                    s=df_pa['parsize'], edgecolors='none', facecolors=df_pa['color'], label='q < 38 au', rasterized=True)
        ax2.scatter(df_pa['a'], df_pa['inc'], alpha=df_pa['alpha'],
                    s=df_pa['parsize'], edgecolors='none', facecolors=df_pa['color'], rasterized=True)
        ax1.text(53.5, 52.5, "Non-physical", fontsize=20,
                 color=opk.BLACK, alpha=0.3, rotation=23, zorder=11, ha='center')
        # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color=opk.BLACK,alpha=0.5,rotation=24)

        ax3.hist(df_highq['a'], bins=nbin, range=(50, 100), density=True,
                 color=opk.LIGHT_BLUE, histtype='step', label='q > 38 au', zorder=2)
        ax3.hist(df_pa['a'][df_pa['q'] < qcut], bins=nbin, range=(50, 100), density=True,
                 color=opk.GRAY, histtype='step', alpha=1, label='q < 38 au', zorder=-10)

        ax3.legend(bbox_to_anchor=(0.97, 0.35), fontsize=14)

        fig.subplots_adjust(hspace=0, right=0.88, left=0.1)
        pos1 = ax1.get_position()
        # print(pos1)
        rect_histy = [pos1.x1, pos1.y0, 0.09, pos1.y1 - pos1.y0]
        ax4 = fig.add_axes(rect_histy, sharey=ax1)

        ax4.hist(df_arange['q'][df_arange['q'] < 38], bins=20, range=(q_low_bound, 38), density=True, orientation='horizontal',
                 color=opk.GRAY, histtype='step', label='q > 38 au', zorder=-1)
        ax4.hist(df_arange['q'][df_arange['q'] > 38], bins=80, range=(38, q_up_bound), density=True, orientation='horizontal',
                 color=opk.LIGHT_BLUE, histtype='step', label='q > 38 au', zorder=-1)
        ax4.tick_params(axis="y", labelleft=False)

        pos1 = ax2.get_position()
        rect_histy = [pos1.x1, pos1.y0, 0.09, pos1.y1 - pos1.y0]
        ax5 = fig.add_axes(rect_histy, sharey=ax2)
        ax5.hist(df_arange['inc'][df_arange['q'] < 38], bins=30, density=True, orientation='horizontal',
                 color=opk.GRAY, histtype='step', label='q < 38 au', zorder=-1)
        ax5.hist(df_arange['inc'][df_arange['q'] > 38], bins=30, density=True, orientation='horizontal',
                 color=opk.LIGHT_BLUE, histtype='step', label='q > 38 au', zorder=-1)
        ax5.tick_params(axis="y", labelleft=False)

        aN = df_pl['a'].iloc[-2]
        if plotREAL:

            kat_data = pd.read_csv("classified-TNOs_q>38_a>49_saved.csv")
            kat_data['q'] = kat_data['a'] * (1-kat_data['e'])
            kat_data = kat_data[kat_data['q'] > 38]

            kat_data['a'] = kat_data['a']/aN_bary*aN
            kat = kat_data[kat_data['a'] > 50].copy()
            kat['marker'] = '^'
            kat.loc[kat['class'] == 'Nresonant', 'marker'] = 'x'
            # print(kat)
            ax1.scatter(kat['a'][kat['marker'] == '^'], kat['q'][kat['marker'] == '^'], alpha=1, s=30,
                        facecolor='none', edgecolor=opk.BLACK, marker='^', lw=1.5, zorder=2)
            ax2.scatter(kat['a'][kat['marker'] == '^'], kat['i'][kat['marker'] == '^'], alpha=1, s=30,
                        facecolor='none', edgecolor=opk.BLACK, marker='^', lw=1.5, zorder=2)
            ax1.scatter(kat['a'][kat['marker'] == 'x'], kat['q'][kat['marker'] == 'x'], alpha=1, s=30,
                        facecolor=opk.BLACK, marker='x', zorder=2)
            ax2.scatter(kat['a'][kat['marker'] == 'x'], kat['i'][kat['marker'] == 'x'], alpha=1, s=30,
                        facecolor=opk.BLACK, marker='x', zorder=2)
            # print(kat['a'], kat['q'])

        markersize = 18
        blue_dot = mlines.Line2D([], [], color=opk.BLUE, marker='.', linestyle='None',
                                 markersize=markersize, label='Detached')
        orange_dot = mlines.Line2D([], [], color=opk.ORANGE, marker='.', linestyle='None',
                                   markersize=markersize, label='Resonant')
        red_dot = mlines.Line2D([], [], color=opk.DARK_RED, marker='.', linestyle='None',
                                markersize=markersize, label='Scattering')
        gray_dot = mlines.Line2D([], [], color=opk.GRAY, marker='.', linestyle='None',
                                 markersize=markersize*0.5, label='q < 38 au')
        black_tri = mlines.Line2D([], [], marker='^', linestyle='None',
                                  markersize=markersize/2.5, label='Real detached', markerfacecolor='none', markeredgecolor=opk.BLACK, markeredgewidth=1.5)
        black_x = mlines.Line2D([], [], marker='X', linestyle='None',
                                markersize=markersize/2, label='Real resonant', color=opk.BLACK, markeredgewidth=0)
        if plotREAL:
            markerlist = [blue_dot, orange_dot,
                          red_dot, gray_dot, black_tri, black_x]
            ax1.legend(loc='upper left', fontsize=14,
                       handles=markerlist, ncol=3).set_zorder(100)
        else:
            markerlist = [blue_dot, orange_dot, red_dot, gray_dot]
            ax1.legend(loc='upper left', fontsize=14,
                       handles=markerlist, ncol=2).set_zorder(100)

        x = np.linspace(inner, outer, 1000)
        y = np.linspace(inner, outer*10000, 1000)
        # ax1.plot(x, x, color='grey', linestyle='dashed', lw=1)
        ax1.fill_between(x, x, y, facecolor=opk.GRAY, zorder=10)

        ax2.set_xlabel("a (au)")
        ax1.set_ylabel("q (au)")
        ax2.set_ylabel("I (deg)")
        ax3.set_ylabel("Fraction")
        ax5.set_xlabel("Fraction")
        fig.align_ylabels([ax1, ax2, ax3])

        ax1.set_xlim(inner, outer)
        ax1.set_ylim(q_low_bound, q_up_bound)

        ax2.set_xlim(inner, outer)
        ax2.set_ylim(i0, i1)

        ax3.set_xlim(inner, outer)
        ax3.set_ylim(0, per1)

        ax4.set_xlim(0, 0.5)
        ax5.set_xlim(0, 0.05)

        ax1.set(xscale='log', yscale='log')

        xticks = np.arange(round(inner/10)*10, outer, 10)
        xticks = np.concatenate(
            (xticks, np.arange(100, round(outer/10)*10+1, 100)))
        # xticks = [50, 60, 70, 80,
        #           90, 100, 200, 300, 400, 500, 600, 700, 800]
        yticks = np.arange(35, round(q_up_bound/5-1)*5+1, 5)
        ax1.set_xticks(xticks)
        ax1.set_yticks(yticks)
        ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        yvals = ax3.get_yticks()
        ax3.set_yticklabels(["{0:.0%}".format(y) for y in yvals])

        ax4.set_xticks([0, 0.25, 0.5])
        ax4.set_xticklabels(["", "25%", "50%"], fontsize=14)
        ax4.xaxis.tick_top()

        ax5.set_xticks([0, 0.025, 0.05])
        ax5.set_xticklabels(["", "2.5%", "5%"], fontsize=14)

        a_N = aN
        klist, jlist, a_rlist = opk.genOuterRes(
            a_N, 49.5, 101, high1=4, high2=50, order_lim=14)
        for k, j, a_res in zip(klist, jlist, a_rlist):
            label = "{0:0d}/{1:0d}".format(j, k)
            ax3.text(a_res, per1+0.015, label,
                     fontsize=k**(-0.4)*14, rotation=60, ha='center')
            # print(a)
            ax1.plot([a_res, a_res], [38, 300], ls='dashed',
                     color=opk.GRAY, zorder=-1, lw=1)
            ax3.plot([a_res, a_res], [0, 1], ls='dashed',
                     color=opk.GRAY, zorder=-1, lw=1)

        hlines = [qcut]
        for line in hlines:
            ax1.plot([inner, outer], [line, line], ls='dashed',
                     color=opk.DARK_RED, zorder=2, alpha=0.8)
            ax4.plot([0, 1], [line, line], ls='solid',
                     color='white', zorder=2, alpha=1)
            ax4.plot([0, 1], [line, line], ls='dashed',
                     color=opk.DARK_RED, zorder=2, alpha=0.8)

        # ax3.set_title("Time: {time:6.2f} Myr".format(
        #     time=t/365.25/1e6), fontsize=20, va='top', y=1.60)

        # fig.tight_layout()
        # fig.subplots_adjust(hspace=0)
        if plotREAL:
            plt.savefig(output_folder +
                        "frame_REAL_{idx:04d}.jpg".format(idx=idx), dpi=200)
            plt.savefig(output_folder +
                        "frame_REAL_{idx:04d}.pdf".format(idx=idx), dpi=200)
        else:
            plt.savefig(output_folder +
                        "frame_{idx:04d}.jpg".format(idx=idx), dpi=200)
            plt.savefig(output_folder +
                        "frame_{idx:04d}.pdf".format(idx=idx), dpi=200)
        print("Saved! frame: {idx:04d}".format(idx=idx))
        plt.close()
        idx += 1


# plot(np.arange(0, 36530000001, 10000000), False)
plot([36530000000], True)
