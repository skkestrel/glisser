import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
import struct
import os
import sys
import pandas as pd
import matplotlib
import orbitplotkit as opk
import random

font = {'size': 18}
matplotlib.rc('font', **font)


# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
# folder1 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNR-Resonant-filter/"
# folder1 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUN-Resonant-filter/"
folder1 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Giant-Synthetic/"
folder2 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Synthetic-4Gyr-filter/"
output_folder = folder1 + "/pics/snapshots_RESO"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
aN_bary = 30.068010875914965
# pd.set_option('display.max_colwidth', None)
# pd.options.display.max_colwidth = 200


def plot(time, noCut, showSimulated, showReal):
    # df_temp = pd.read_csv(
    #     folder1+"classes.txt", names=['id', 'RESO', 'k1', 'k2'], delimiter=' ').set_index('id')

    # time = [0, 21520000000]
    idx = 1
    # time = [0]

    par_size = 3
    qcut = 38
    q_low_bound = 37
    q_up_bound = 64
    i0, i1 = 0, 62
    per1 = 1.10
    a0, a1 = 50, 100
    inner, outer = a0-0.5, a1+1

    nbin = 200
    lw = 2
    s = 40

    for t in time:
        fig, [ax3, ax1, ax2] = plt.subplots(3, figsize=(
            11, 7.5), gridspec_kw={'height_ratios': [1, 4, 1]}, sharex=True)

        pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(0)
        df_pl = pd.read_csv(
            folder2+pl_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')
        df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
        df_pl['inc'] = np.rad2deg(df_pl['inc'])
        df_pl['Omega'] = np.rad2deg(df_pl['Omega'])
        df_pl['omega'] = np.rad2deg(df_pl['omega'])
        df_pl['M'] = np.rad2deg(df_pl['M'])


        if not noCut:
            pa_txt = "_rogue_sim_export_4Gyr_q35_detached.csv"
        else:
            pa_txt = "_rogue_sim_export_4Gyr_detached.csv"
        df_pa = pd.read_csv(pa_txt, delimiter=',').set_index('id')
        # df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])

        # pa_txt = "reoutput/snapshots_RESO/temp_{0:d}.csv".format(0)
        # df_pa0 = pd.read_csv(folder1+pa_txt, delimiter=',').set_index('id')
        # df_pa0['q'] = df_pa0['a'] * (1 - df_pa0['e'])
        df_pa['inc'] = np.rad2deg(df_pa['inc'])
        # df_pa['Omega'] = np.rad2deg(df_pa['Omega'])
        # df_pa['omega'] = np.rad2deg(df_pa['omega'])
        # df_pa['M'] = np.rad2deg(df_pa['M'])
        # df_pa['q0'] = df_pa0['q']

        # if not noCut:
        #     df_pa = df_pa[df_pa['q0'] < 35]

        # print(df_pa)

        df_pa['alpha'] = 1
        # df_pa.loc[df_pa['RESO'] == 1, 'alpha'] = 0.5
        df_pa.loc[df_pa['q'] < qcut, 'alpha'] = 0
        # df_pa.loc[df_pa['q'] > 40, 'alpha'] = 1

        df_pa['parsize'] = 1.5
        df_pa.loc[df_pa['q'] > qcut, 'parsize'] = 5
        # df_pa.loc[df_pa['q'] > 40, 'parsize'] = 4

        df_pa['color'] = opk.GRAY

        df_highq = df_pa[df_pa['q'] > qcut].copy()
        df_arange = df_highq[df_highq['a'].between(a0, a1)].copy()
        # df_used = df_arange[df_arange['RESO'] == 0]
        df_used = df_arange
        # df_highq = df_highq[df_highq['a'].between(a0, a1)]
        # print(df_highq)
        print(df_used)

        # df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
        # ax1.scatter(df_pl['a'], df_pl['q'], alpha=0.8, s=100,
        #             color=opk.DARK_RED, marker='x', lw=2, zorder=2)
        # ax2.scatter(df_pl['a'], np.rad2deg(df_pl['inc']), alpha=0.8,
        #             s=100, color=opk.DARK_RED, marker='x', lw=2, zorder=2)

        ax1.scatter(df_used['a'], df_used['q'], alpha=df_used['alpha'], rasterized=True,
                    s=df_used['parsize'], edgecolors='none', facecolors=df_used['color'], label='q < 38 au')
        ax2.scatter(df_used['a'], df_used['inc'], alpha=df_used['alpha'], rasterized=True,
                    s=df_used['parsize'], edgecolors='none', facecolors=df_used['color'])
        ax1.text(54.5, 52.5, "Non-physical", fontsize=18,
                 color=opk.BLACK, alpha=0.3, rotation=36, zorder=11, ha='center')
        # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color=opk.BLACK,alpha=0.5,rotation=24)

        # ax3.hist(df_pa['a'][df_pa['q'] < qcut], bins=nbin, range=(a0, a1), density=True,
        #          color=opk.GRAY, histtype='step', alpha=1, label='q < 38 au', zorder=-10)

        fig.subplots_adjust(hspace=0, right=0.83, left=0.15)
        pos1 = ax1.get_position()
        # print(pos1)
        rect_histy = [pos1.x1, pos1.y0, 0.09, pos1.y1 - pos1.y0]
        ax4 = fig.add_axes(rect_histy, sharey=ax1)

        pos1 = ax2.get_position()
        rect_histy = [pos1.x1, pos1.y0, 0.09, pos1.y1 - pos1.y0]
        ax5 = fig.add_axes(rect_histy, sharey=ax2)

        aN = df_pl['a'].iloc[-1]
        if showReal:
            kat_data = pd.read_csv("classified-TNOs_q>38_a>49_saved.csv")
            kat_data['q'] = kat_data['a'] * (1-kat_data['e'])
            kat_data = kat_data[kat_data['q'] > 38]

            kat_data['a'] = kat_data['a']/aN_bary*aN
            kat = kat_data[kat_data['a'].between(a0, a1)].copy()
            kat['marker'] = '^'
            kat.loc[kat['class'] == 'Nresonant', 'marker'] = 'x'
            kat.loc[kat['class'] == 'scattering', 'marker'] = 'D'
            kat_detached = kat[kat['marker'] == '^']
            kat_resonant = kat[kat['marker'] == 'x']
            num_real = kat.shape[0]
            num_det = kat_detached.shape[0]
            num_res = kat_resonant.shape[0]
            print(kat_detached.sort_values(by=['i']))
            print(num_det, num_res)
            print(num_det / num_real, num_res/num_real)

            ax1.scatter(kat_detached['a'], kat_detached['q'], alpha=1, s=s,
                        facecolor='none', edgecolor=opk.BLACK, marker='^', lw=lw, zorder=2)
            ax2.scatter(kat_detached['a'], kat_detached['i'], alpha=1, s=s,
                        facecolor='none', edgecolor=opk.BLACK, marker='^', lw=lw, zorder=2)
            ax3.hist(kat_detached['a'], bins=nbin, range=(a0, a1), density=True,  cumulative=True,
                     color=opk.BLACK, histtype='step', label='Kat', zorder=2, lw=2)
            ax4.hist(kat_detached['q'], bins=nbin, range=(qcut, q_up_bound), density=True, orientation='horizontal',  cumulative=True,
                     color=opk.BLACK, histtype='step', label='Kat', zorder=-1, lw=2)
            ax5.hist(kat_detached['i'], bins=nbin, range=(i0, i1), density=True, orientation='horizontal', cumulative=True,
                     color=opk.BLACK, histtype='step', label='Kat', zorder=-1, lw=2)

            # real_data = pd.read_csv("all_distant_tnos_classified.csv")
            # real_data = real_data[real_data['Source'] == 'KAT']
            # real_data = real_data[real_data['q'] > 38]
            # real_data['a'] = real_data['a']/aN_bary*aN
            # real = real_data[real_data['a'].between(a0, a1)].copy()
            # real_detached = real[real['RESO'] == 0]
            # real_resonant = real[real['RESO'] > 0]
            # num_real = real.shape[0]
            # num_det = real_detached.shape[0]
            # num_res = real_resonant.shape[0]
            # print(real_detached)
            # print(num_det, num_res)
            # print(num_det / num_real, num_res/num_real)

            # ax1.scatter(real_detached['a'], real_detached['q'], alpha=1, s=s,
            #             facecolor='none', edgecolor=opk.BLACK, marker='^', lw=lw, zorder=2)
            # ax2.scatter(real_detached['a'], real_detached['i'], alpha=1, s=s,
            #             facecolor='none', edgecolor=opk.BLACK, marker='^', lw=lw, zorder=2)
            # ax3.hist(real_detached['a'], bins=nbin, range=(a0, a1), density=True,  cumulative=True,
            #          color=opk.BLACK, histtype='step', label='Kat', zorder=2, lw=2)
            # ax4.hist(real_detached['q'], bins=nbin, range=(qcut, q_up_bound), density=True, orientation='horizontal',  cumulative=True,
            #          color=opk.BLACK, histtype='step', label='Kat', zorder=-1, lw=2)
            # ax5.hist(real_detached['i'], bins=nbin, range=(i0, i1), density=True, orientation='horizontal', cumulative=True,
            #          color=opk.BLACK, histtype='step', label='Kat', zorder=-1, lw=2)

        if showSimulated:
            if noCut:
                readfile = "q0nocut_detections_4Gyr_detached.txt"
            else:
                readfile = "q0le35_detections_4Gyr_detached.txt"
            matthew_data = pd.read_csv(readfile, names=[
                'a', 'e', 'i', 'q', 'r', 'M', 'node', 'peri', 'm_rand', 'H_rand', 'color', 'Comments'], skiprows=2, header=None, delim_whitespace=True)
            matthew_data = matthew_data[matthew_data['q'] > qcut]
            # matthew_data['a'] = matthew_data['a']/aN_bary*aN
            matthew = matthew_data[matthew_data['a'].between(a0, a1)].copy()
            matthew = matthew
            matthew_sample = matthew.sample(num_det)
            print(matthew)

            ax1.scatter(matthew_sample['a'], matthew_sample['q'], alpha=1, s=s,
                        color=opk.DARK_RED, marker='x', lw=lw, zorder=1)
            ax2.scatter(matthew_sample['a'], matthew_sample['i'], alpha=1, s=s,
                        color=opk.DARK_RED, marker='x', lw=lw, zorder=1)
            ax3.hist(matthew['a'], bins=nbin, range=(a0, a1), density=True,  cumulative=True,
                     color=opk.DARK_RED, histtype='step', label='Matthew', zorder=2, lw=2)
            ax4.hist(matthew['q'], bins=nbin, range=(qcut, q_up_bound), density=True, orientation='horizontal',  cumulative=True,
                     color=opk.DARK_RED, histtype='step', label='Matthew', zorder=-1, lw=2)
            ax5.hist(matthew['i'], bins=nbin, range=(i0, i1), density=True, orientation='horizontal', cumulative=True,
                     color=opk.DARK_RED, histtype='step', label='Matthew', zorder=-1, lw=2)

        ax3.hist(df_used['a'], bins=nbin, range=(a0, a1), density=True,  cumulative=True,
                 color=opk.GRAY, histtype='step', label='Simulated', zorder=0, lw=2)
        ax4.hist(df_used['q'], bins=nbin, range=(qcut, q_up_bound), density=True, orientation='horizontal',  cumulative=True,
                 color=opk.GRAY, histtype='step', label='Simulated', zorder=-5, lw=2)
        ax4.tick_params(axis="y", labelleft=False)
        ax5.hist(df_used['inc'], bins=nbin, range=(i0, i1), density=True, orientation='horizontal', cumulative=True,
                 color=opk.GRAY, histtype='step', label='Simulated', zorder=-5, lw=2)
        ax5.tick_params(axis="y", labelleft=False)

        markersize = 18
        blue_dot = mlines.Line2D([], [], color=opk.BLUE, marker='.', linestyle='None',
                                 markersize=markersize, label='Detached')
        orange_dot = mlines.Line2D([], [], color=opk.ORANGE, marker='.', linestyle='None',
                                   markersize=markersize, label='Resonant')
        red_dot = mlines.Line2D([], [], color=opk.DARK_RED, marker='.', linestyle='None',
                                markersize=markersize, label='Scattering')

        if noCut:
            label = 'Intrinsic detached'
        else:
            label = 'Intrinsic detached'
        gray_dot = mlines.Line2D([], [], color=opk.GRAY, marker='.', linestyle='None',
                                 markersize=markersize*0.5, label=label)
        black_tri = mlines.Line2D([], [], marker='^', linestyle='None',
                                  markersize=markersize/2.5, label='Real detached', markerfacecolor='none', markeredgecolor=opk.BLACK, markeredgewidth=1.5)
        red_x = mlines.Line2D([], [], marker='X', linestyle='None',
                              markersize=markersize/2, label='Simulated detections', color=opk.DARK_RED, markeredgewidth=0)

        red_line = matplotlib.lines.Line2D(
            [], [], c=opk.RED, label='Simulated')
        black_line = matplotlib.lines.Line2D([], [], c=opk.BLACK, label='Real')
        gray_line = matplotlib.lines.Line2D(
            [], [], c=opk.GRAY, label='Intrinsic')

        markerlist = [gray_dot]
        linelist = [gray_line]

        if showSimulated:
            markerlist.append(red_x)
            linelist.append(red_line)

        if showReal:
            markerlist.append(black_tri)
            linelist.append(black_line)
        ax1.legend(loc='upper center', fontsize=13,
                   handles=markerlist, ncol=1).set_zorder(100)
        ax3.legend(bbox_to_anchor=(1, 1.30), fontsize=13,
                   handles=linelist).set_zorder(100)

        x = np.linspace(inner, outer, 1000)
        y = np.linspace(inner, outer*10000, 1000)
        ax1.fill_between(x, x, y, facecolor=opk.GRAY, zorder=10)

        ax2.set_xlabel("a (au)")
        ax1.set_ylabel("q (au)")
        ax2.set_ylabel("I (deg)")
        ax3.set_ylabel("Cumulative\nfraction")
        ax5.set_xlabel("Cumulative\nfraction")
        fig.align_ylabels([ax1, ax2, ax3])

        ax1.set_xlim(inner, outer)
        ax1.set_ylim(q_low_bound, q_up_bound)

        ax2.set_xlim(inner, outer)
        ax2.set_ylim(i0-3, i1)

        ax3.set_xlim(inner, outer)
        ax3.set_ylim(0, per1)

        ax4.set_xlim(0, per1)
        ax5.set_xlim(0, per1)

        ax1.set(xscale='log', yscale='log')

        xticks = np.arange(round(inner/10)*10, outer, 10)
        xticks = np.concatenate(
            (xticks, np.arange(100, round(outer/10)*10+1, 100)))
        yticks = np.arange(40, round(q_up_bound/5-1)*5+1, 5)
        yticks = np.append(yticks, 38)
        ax1.set_xticks(xticks)
        ax1.set_yticks(yticks)
        ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        yvals = ax3.get_yticks()
        ax3.set_yticklabels(["{0:.0%}".format(y) for y in yvals])

        ax4.set_xticks([0, 1.0])
        ax4.set_xticklabels(["", "100%"], fontsize=14)
        ax4.xaxis.tick_top()

        ax5.set_xticks([0, 1.0])
        ax5.set_xticklabels(["", "100%"], fontsize=14)

        a_N = aN
        klist, jlist, a_rlist = opk.genOuterRes(
            a_N, inner, outer, high1=3, high2=50, order_lim=8)
        for k, j, a_res in zip(klist, jlist, a_rlist):
            label = "{0:0d}/{1:0d}".format(j, k)
            ax3.text(a_res, per1+0.15, label,
                     fontsize=k**(-0.4)*14, rotation=60, ha='center')
            # print(a)
            ax1.plot([a_res, a_res], [38, 300], ls='dashed',
                     color=opk.GRAY, zorder=-1, lw=1)
            ax3.plot([a_res, a_res], [0, per1], ls='dashed',
                     color=opk.GRAY, zorder=-1, lw=1)

        hlines = [qcut]
        for line in hlines:
            ax1.plot([inner, outer], [line, line], ls='dashed',
                     color=opk.DARK_RED, zorder=2, alpha=0.8)
            ax4.plot([0, per1], [line, line], ls='solid',
                     color='white', zorder=2, alpha=1)
            ax4.plot([0, per1], [line, line], ls='dashed',
                     color=opk.DARK_RED, zorder=2, alpha=0.8)

        # ax3.set_title("Time: {time:6.2f} Myr".format(
        #     time=t/365.25/1e6), fontsize=20, va='top', y=1.60)

        # fig.tight_layout()
        # fig.subplots_adjust(hspace=0)
        if noCut:
            string = 'nocut'
        else:
            string = 'q<35'
        name = '_detection_'+string+'_4Gyr_ref0'
        if showReal:
            name = '_detection_'+string+'_4Gyr_ref1'
        if showSimulated:
            name = '_detection_'+string+'_4Gyr_all'
        plt.savefig(
            name+'.jpg', dpi=200)
        plt.savefig(
            name+'.pdf', dpi=200)
        print("Saved! frame: {name}.jpg".format(name=name))
        print("Saved! frame: {name}.pdf".format(name=name))
        plt.close()
        idx += 1


time = [35020000000]

# plot(time, False, False, False)
# plot(time, False, False, True)
# plot(time, True, True, True)
plot(time, False, True, True)
# plot(time, False, False, True)
