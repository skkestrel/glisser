import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import struct
import os
import pandas as pd
import matplotlib
import orbitplotkit as opk

# font = {'weight' : 'bold',
#         'size'   : 16}
# matplotlib.rc('font', **font)


# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
# folder1 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNR-Resonant-filter/"
# folder1 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUN-Resonant-filter/"
folder1 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Synthetic-1Gyr-filter/"
output_folder = "/home/yhuang/GLISSER/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# df_temp = pd.read_csv(
#     folder1+"classes.txt", names=['id', 'RESO', 'k1', 'k2'], delimiter=' ').set_index('id')
time = np.arange(0, 36500000001, 50000000)
time = [0]
idx = 1
# time = [0]

par_size = 1


for t in time:
    pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
    df_pl = pd.read_csv(
        folder1+pl_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')
    df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
    df_pl['inc'] = np.rad2deg(df_pl['inc'])
    df_pl['Omega'] = np.rad2deg(df_pl['Omega'])
    df_pl['omega'] = np.rad2deg(df_pl['omega'])
    df_pl['M'] = np.rad2deg(df_pl['M'])

    pa_txt = "reoutput/snapshots/temp_{0:d}.csv".format(t)
    df_pa = pd.read_csv(folder1+pa_txt, delimiter=',').set_index('id')
    df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])
    # df_pa['q0'] = df_pa0['q']
    # df_pa['inc0'] = df_pa0['inc']
    df_pa['inc'] = np.rad2deg(df_pa['inc'])
    df_pa['Omega'] = np.rad2deg(df_pa['Omega'])
    df_pa['omega'] = np.rad2deg(df_pa['omega'])
    df_pa['M'] = np.rad2deg(df_pa['M'])
    print(df_pa)

    fig, [ax1, ax2] = plt.subplots(2, figsize=(
        15, 8), gridspec_kw={'height_ratios': [4, 1]}, sharex=True)

    qcut = 38

    df_pa['alpha'] = 1
    df_pa.loc[df_pa['RESO'] == 1, 'alpha'] = 0.5
    df_pa.loc[df_pa['q'] < qcut, 'alpha'] = 0.35
    # df_pa.loc[df_pa['q'] > 40, 'alpha'] = 1

    df_pa['parsize'] = 1
    df_pa.loc[df_pa['q'] > qcut, 'parsize'] = 1
    # df_pa.loc[df_pa['q'] > 40, 'parsize'] = 4

    df_pa['color'] = 'gray'
    df_pa.loc[df_pa['RESO'] == 0, 'color'] = opk.LIGHT_BLUE
    df_pa.loc[df_pa['q'] > qcut, 'color'] = opk.BLUE
    df_pa.loc[df_pa['RESO'] == -1, 'color'] = opk.RED
    df_pa.loc[df_pa['RESO'] > 1, 'color'] = opk.ORANGE

    df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
    ax1.scatter(df_pl['a'], df_pl['q'], alpha=0.8, s=100,
                color=opk.RED, marker='x', lw=2, zorder=2)
    ax2.scatter(df_pl['a'], np.rad2deg(df_pl['inc']), alpha=0.8,
                s=100, color=opk.RED, marker='x', lw=2, zorder=2)

    ax1.scatter(df_pa['a'], df_pa['q'], alpha=df_pa['alpha'],
                s=df_pa['parsize'], edgecolors='none', facecolors=df_pa['color'])
    ax2.scatter(df_pa['a'], df_pa['inc'], alpha=df_pa['alpha'],
                s=df_pa['parsize'], edgecolors='none', facecolors=df_pa['color'])
    ax1.text(53, 58, "Non-physical", fontsize=18,
             color='black', alpha=0.25, rotation=38, zorder=11)
    # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

    aN = df_pl['a'].iloc[-2]
    outer = 800
    inner = 49
    x = np.linspace(inner, outer, 1000)
    y = np.linspace(inner, outer*10000, 1000)
    # ax1.plot(x, x, color='grey', linestyle='dashed', lw=1)
    ax1.fill_between(x, x, y, facecolor='#BEBEBE', zorder=10)

    # ax1.plot([30, outer], [38, 38], color=red, alpha=0.5, linestyle='dashed',lw=2)
    # ax2.plot([30, outer], [6, 6], color=red, linestyle='dashed',lw=2)
    # ax1.plot([80, outer], [80, 80], color=red, alpha=0.3, linestyle='dashed',lw=2)
    # for pos,label in zip([aN,38], ["q = 30", "q = 38"]):
    #     ax1.text(15, pos-2, label, fontsize = 12, color=red)

    # rect2 = patches.Rectangle((4, 80), 46, 20, edgecolor='none', facecolor='green', alpha=0.2)
    # ax2.add_patch(rect2)

    q_low_bound = 29.5
    q_up_bound = 150

    ax2.set_xlabel("a (au)")
    ax1.set_ylabel("q (au)")
    ax2.set_ylabel("I (deg)")

    ax1.set_xlim(inner, outer)
    ax1.set_ylim(q_low_bound, q_up_bound)
    ax2.set_xlim(inner, outer)
    ax2.set_ylim(0, 62)

    ax1.set(xscale='log', yscale='log')

    # xticks = np.arange(round(inner/100)*100, outer, 100)
    # xticks = np.concatenate(
    #     (xticks, np.arange(100, round(outer/100)*100+1, 100)))
    xticks = [50, 60, 70, 80,
              90, 100, 200, 300, 400, 500, 600]
    yticks = np.arange(round(q_low_bound/10)*10, round(q_up_bound/10)*10+1, 10)
    ax1.set_xticks(xticks)
    ax1.set_yticks(yticks)
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    a_N = aN

    kat_data = pd.read_csv("classified-TNOs_q>38_a>49_saved.csv")
    kat_data['q'] = kat_data['a'] * (1-kat_data['e'])
    kat_data = kat_data[kat_data['q'] > 38]

    kat = kat_data[kat_data['a'] > 50].copy()
    kat['marker'] = '^'
    kat.loc[kat['class'] == 'Nresonant', 'marker'] = 'x'
    kat.loc[kat['class'] == 'scattering', 'marker'] = 'D'
    # print(kat)
    ax1.scatter(kat['a'][kat['marker'] == '^'], kat['q'][kat['marker'] == '^'], alpha=1, s=30,
                facecolor='none', edgecolor=opk.BLACK, marker='^', lw=1.5, zorder=2)
    ax2.scatter(kat['a'][kat['marker'] == '^'], kat['i'][kat['marker'] == '^'], alpha=1, s=30,
                facecolor='none', edgecolor=opk.BLACK, marker='^', lw=1.5, zorder=2)
    ax1.scatter(kat['a'][kat['marker'] == 'x'], kat['q'][kat['marker'] == 'x'], alpha=1, s=30,
                facecolor=opk.BLACK, marker='x', zorder=2)
    ax2.scatter(kat['a'][kat['marker'] == 'x'], kat['i'][kat['marker'] == 'x'], alpha=1, s=30,
                facecolor=opk.BLACK, marker='x', zorder=2)
    ax1.scatter(kat['a'][kat['marker'] == 'D'], kat['q'][kat['marker'] == 'D'], alpha=1, s=30,
                facecolor=opk.BLACK, marker='D', zorder=2)
    ax2.scatter(kat['a'][kat['marker'] == 'D'], kat['i'][kat['marker'] == 'D'], alpha=1, s=30,
                facecolor=opk.BLACK, marker='D', zorder=2)
    klist, jlist, a_rlist = opk.genOuterRes(
        a_N, 48, 120, high1=2, high2=100, order_lim=25)
    for k, j, a_res in zip(klist, jlist, a_rlist):
        label = "{0}/{1}".format(j, k)
        ax1.text(a_res-0.8, q_up_bound+1, label, fontsize=9, rotation=70)
        # print(a)
        ax1.plot([a_res, a_res], [38, 600], ls='dashed',
                 color='#BEBEBE', zorder=-1, lw=1)

    hlines = [qcut]
    for line in hlines:
        ax1.plot([inner, outer], [line, line], ls='dashed',
                 color=opk.RED, zorder=2, alpha=0.8)

    plt.title("Time: {time:6.3f} Myr".format(time=1000))
    fig.subplots_adjust(hspace=0)
    plt.savefig(output_folder +
                "z_figure_{idx:04d}.jpg".format(idx=idx), dpi=250)
    print("Saved! frame: {idx:04d}".format(idx=idx))
    plt.close()
    idx += 1
