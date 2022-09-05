import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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


# jpl_data = pd.read_csv("sbdb_query_results_1.19.csv")
# jpl_data = jpl_data[jpl_data['q'] > 38]
# jpl = jpl_data[jpl_data['a'] > 50].copy()
# print(jpl)

# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Full-Scattering/"
output_folder = folder + "pics/snapshots/zoomout/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# time = [0, 182000000000]
# time = [6275000000]
time = np.arange(0, 36525000001, 25000000)
# print(time.shape[0])
idx = 1
time = [0, 10850000000]
par_size = 8

for t in time:
    pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
    df_pl = pd.read_csv(
        folder+pl_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ')
    # print(df_pl)

    pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
    df_pa = pd.read_csv(
        folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ')
    df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
    df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])

    # df_pa2 = df_pa[df_pa['id'] > 90000]
    # df_pa = df_pa[df_pa['id'] <= 90000]

    fig, [ax1, ax2] = plt.subplots(2, figsize=(
        18, 10), gridspec_kw={'height_ratios': [3, 1]})

    # ax1.scatter(jpl['a'], jpl['q'], alpha=0.5, s=30,
    #             facecolor='none', edgecolor=red, marker='^', lw=2, zorder=2)
    # ax2.scatter(jpl['a'], jpl['i'], alpha=0.5, s=30,
    #             facecolor='none', edgecolor=red, marker='^', lw=2, zorder=2)

    ax1.scatter(df_pl['a'], df_pl['q'], alpha=0.8, s=150,
                color='black', marker='x', lw=3, zorder=2)
    ax2.scatter(df_pl['a'], np.rad2deg(df_pl['inc']), alpha=0.8,
                s=150, color='black', marker='x', lw=3, zorder=2)

    ax1.scatter(df_pa['a'][df_pa['q'] > 38], df_pa['q'][df_pa['q'] > 38],
                alpha=0.8, s=par_size, edgecolors='none', facecolors='C0')
    ax1.scatter(df_pa['a'][df_pa['q'] < 38], df_pa['q'][df_pa['q'] < 38],
                alpha=0.15, s=par_size, edgecolors='none', facecolors='C0')

    ax2.scatter(df_pa['a'][df_pa['q'] > 38], np.rad2deg(
        df_pa['inc'][df_pa['q'] > 38]), alpha=0.8, s=par_size, edgecolors='none', facecolors='C0')
    ax2.scatter(df_pa['a'][df_pa['q'] < 38], np.rad2deg(
        df_pa['inc'][df_pa['q'] < 38]), alpha=0.1, s=par_size, edgecolors='none', facecolors='C0')
    # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
    # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

    aN = df_pl['a'].iloc[3]
    outer = 500
    inner = 28
    x = np.linspace(inner, outer, 1000)
    y = np.linspace(inner, outer*10000, 1000)
    ax1.plot(x, x, color='grey', linestyle='dashed', lw=1)
    ax1.fill_between(x, x, y, facecolor='gray', alpha=0.5)

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

    ax2.set_xlabel("a (au)")
    ax1.set_ylabel("q (au)")
    ax2.set_ylabel("I (deg)")

    # ax1.set_xscale("log")
    # ax1.set_yscale("log")
    # ax2.set_xscale("log")
    inner_bound = inner
    outer_bound = outer
    ax1.set_xlim(inner_bound, outer_bound)
    ax2.set_xlim(inner_bound, outer_bound)

    q_low_bound = 25
    q_up_bound = 120
    ax1.set_ylim(q_low_bound, q_up_bound)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    xticks = [30, 40, 50, 60, 70, 80, 90 ,100, 200 ,300 ,400, 500]
    yticks = np.arange(30,120,10)
    ax1.set_xticks(xticks)
    ax1.set_yticks(yticks)
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    ax2.set_ylim(0, 60)
    ax2.set_xscale("log")
    ax2.set_xticks(xticks)
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    a_N = aN
    klist, jlist, a_rlist = opk.genOuterRes(
        a_N, 49, 126, high1=2, high2=50, order_lim=45)
    for k, j, a_res in zip(klist, jlist, a_rlist):
        label = "{0}/{1}".format(j, k)
        ax1.text(a_res-0.3, q_up_bound+1, label, fontsize=8, rotation=70)
        # print(a)
        ax1.plot([a_res, a_res], [0, 300], ls='dashed',
                 color='gray', zorder=1, alpha=0.6)

    hlines = [30, 38]
    for line in hlines:
        ax1.plot([inner, outer], [line, line], ls='dashed',
                 color='red', zorder=2, alpha=0.6)

    plt.title("Time: {time:6.3f} Myr".format(time=t/365.25/1e6))
    plt.tight_layout()

    plt.savefig(output_folder +
                "frame_{idx:04d}_zoomout.jpg".format(idx=idx), dpi=200)
    print("Saved! frame: {idx:04d}".format(idx=idx))
    plt.close()
    idx += 1
