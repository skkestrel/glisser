import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import struct
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


# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNR-Resonant-lowq/"
output_folder = folder + "pics/snapshots/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# time = [0, 182000000000]
# time = [6275000000]
time = np.arange(0, 36500000001, 50000000)
idx = 1
time = [0, 36520000000]


for t in time:
    pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
    df_pl = pd.read_csv(
        folder+pl_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ')
    print(df_pl)

    pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
    df_pa = pd.read_csv(
        folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'If'], delimiter=' ')
    print(df_pa)

    fig, [ax1, ax2] = plt.subplots(2, figsize=(9, 10))
    df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
    ax1.scatter(df_pl['a'], df_pl['q'], alpha=0.8, s=100,
                color=red, marker='x', lw=2, zorder=2)
    ax2.scatter(df_pl['a'], np.rad2deg(df_pl['inc']), alpha=0.8,
                s=100, color=red, marker='x', lw=2, zorder=2)

    df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])
    par_size = 2
    ax1.scatter(df_pa['a'], df_pa['q'], alpha=0.4,
                s=par_size, edgecolors='none', facecolors='C0')
    ax2.scatter(df_pa['a'], np.rad2deg(df_pa['If']), alpha=0.4,
                s=par_size, edgecolors='none', facecolors='C0')
    # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
    # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

    aN = df_pl['a'].iloc[-1]
    outer = 50
    inner = 30
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
    ax2.set_ylabel("If (deg)")

    # ax1.set_xscale("log")
    # ax1.set_yscale("log")
    # ax2.set_xscale("log")
    ax1.set_xlim(40, 48)
    ax1.set_ylim(30, 48)

    ax2.set_xlim(40, 48)
    ax2.set_ylim(0.1, 40)
    ax2.set_yscale("log")

    a_N = aN
    a_32 = opk.resonanceLocationBYA(a_N, 2, 3)
    a_85 = opk.resonanceLocationBYA(a_N, 5, 8)
    a_53 = opk.resonanceLocationBYA(a_N, 3, 5)
    a_127 = opk.resonanceLocationBYA(a_N, 7, 12)
    a_74 = opk.resonanceLocationBYA(a_N, 4, 7)
    a_95 = opk.resonanceLocationBYA(a_N, 5, 9)
    a_21 = opk.resonanceLocationBYA(a_N, 1, 2)

    for a, label in zip([a_32, a_85, a_53, a_127, a_74, a_95, a_21],
                        [" 3/2", " 8/5", " 5/3", "12/7", " 7/4", " 9/5", " 2/1"]):
        ax1.annotate(label, fontsize=12,
                     xy=(a, 48), xycoords='data',
                     xytext=(a-0.1, 49), textcoords='data',
                     arrowprops=dict(arrowstyle="-",
                                     connectionstyle="arc3"),
                     )

    plt.title("Time: {time:6.3f} Myr".format(time=t/365.25/1e6))

    plt.savefig(output_folder+"frame_{idx:04d}.jpg".format(idx=idx), dpi=200)
    print("Saved! frame: {idx:04d}".format(idx=idx))
    plt.close()
    idx += 1
