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


# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNR-Resonant-4Gy-filter/"
output_folder = folder + "pics/snapshots/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

df = pd.read_csv(folder + "classes.txt",
                 names=['id', 'RESO', 'k', 'j'], delimiter=' ').set_index('id')

time = [0]
# time = [6275000000]
# time = np.arange(0, 36520000001, 20000000)
# print(time.shape[0])
idx = 1
# time = [0, 3640000000]
par_size = 8

for t in time:
    pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
    df_pl = pd.read_csv(
        folder+pl_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ').set_index('id')
    # print(df_pl)

    pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
    df_pa = pd.read_csv(
        folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ').set_index('id')
    # df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
    # df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])
    df_pa['M'] = opk.F2M(df_pa['f'], df_pa['e'])
    # del df_pa['f']
    df_pa['RESO'] = df['RESO']
    df_pa['k'] = df['k']
    df_pa['j'] = df['j']
    df_pa.to_csv("JSUNR_4Gyr.csv")

    print(df_pa)

    df_sca = df_pa[df_pa['RESO'] < 0]
    df_det = df_pa[df_pa['RESO'] == 0]
    df_res = df_pa[df_pa['RESO'] > 0]

    # print(df_sca.count(), df_det.count(), df_res.count(), df_pa.count())

    # df_pa2 = df_pa[df_pa['id'] > 90000]
    # df_pa = df_pa[df_pa['id'] <= 90000]

    fig, [ax1, ax2] = plt.subplots(2, figsize=(
        18, 10), gridspec_kw={'height_ratios': [3, 1]})

    # ax1.scatter(df_pl['a'], df_pl['q'], alpha=0.8, s=150,
    #             color='black', marker='x', lw=3, zorder=2)
    # ax2.scatter(df_pl['a'], np.rad2deg(df_pl['inc']), alpha=0.8,
    #             s=150, color='black', marker='x', lw=3, zorder=2)

    ax1.scatter(df_sca['a'][df_sca['q'] > 38], df_sca['q'][df_sca['q'] > 38],
                alpha=0.8, s=par_size, edgecolors='none', facecolors='C1')
    ax1.scatter(df_sca['a'][df_sca['q'] < 38], df_sca['q'][df_sca['q'] < 38],
                alpha=0.15, s=par_size, edgecolors='none', facecolors='C1')

    ax1.scatter(df_det['a'][df_det['q'] > 38], df_det['q'][df_det['q'] > 38],
                alpha=0.8, s=par_size, edgecolors='none', facecolors='C0')
    ax1.scatter(df_det['a'][df_det['q'] < 38], df_det['q'][df_det['q'] < 38],
                alpha=0.15, s=par_size, edgecolors='none', facecolors='C0')

    ax1.scatter(df_res['a'][df_res['q'] > 38], df_res['q'][df_res['q'] > 38],
                alpha=0.8, s=par_size, edgecolors='none', facecolors='C2')
    ax1.scatter(df_res['a'][df_res['q'] < 38], df_res['q'][df_res['q'] < 38],
                alpha=0.15, s=par_size, edgecolors='none', facecolors='C2')

    ax2.scatter(df_sca['a'][df_sca['q'] > 38], np.rad2deg(
        df_sca['inc'][df_sca['q'] > 38]), alpha=0.8, s=par_size, edgecolors='none', facecolors='C1')
    ax2.scatter(df_sca['a'][df_sca['q'] < 38], np.rad2deg(
        df_sca['inc'][df_sca['q'] < 38]), alpha=0.1, s=par_size, edgecolors='none', facecolors='C1')

    ax2.scatter(df_det['a'][df_det['q'] > 38], np.rad2deg(
        df_det['inc'][df_det['q'] > 38]), alpha=0.8, s=par_size, edgecolors='none', facecolors='C0')
    ax2.scatter(df_det['a'][df_det['q'] < 38], np.rad2deg(
        df_det['inc'][df_det['q'] < 38]), alpha=0.1, s=par_size, edgecolors='none', facecolors='C0')

    ax2.scatter(df_res['a'][df_res['q'] > 38], np.rad2deg(
        df_res['inc'][df_res['q'] > 38]), alpha=0.8, s=par_size, edgecolors='none', facecolors='C2')
    ax2.scatter(df_res['a'][df_res['q'] < 38], np.rad2deg(
        df_res['inc'][df_res['q'] < 38]), alpha=0.1, s=par_size, edgecolors='none', facecolors='C2')
    # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
    # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

    aN = df_pl['a'].iloc[3]
    outer = 126
    inner = 49
    x = np.linspace(inner, outer, 1000)
    y = np.linspace(inner, outer*10000, 1000)
    ax1.plot(x, x, color='grey', linestyle='dashed', lw=1)
    ax1.fill_between(x, x, y, facecolor='gray', alpha=0.5)

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

    q_low_bound = 28
    q_up_bound = 75
    ax1.set_ylim(q_low_bound, q_up_bound)
    ax2.set_ylim(0, 60)
    # ax2.set_yscale("log")

    a_N = aN
    klist, jlist, a_rlist = opk.genOuterRes(
        a_N, 49, 126, high1=5, high2=50, order_lim=45)
    for k, j, a_res in zip(klist, jlist, a_rlist):
        label = "{0}/{1}".format(j, k)
        ax1.text(a_res-0.3, q_up_bound+1, label, fontsize=8, rotation=70)
        # print(a)
        ax1.plot([a_res, a_res], [0, 100], ls='dashed',
                 color='gray', zorder=1, alpha=0.6)

    hlines = [30, 38]
    for line in hlines:
        ax1.plot([inner, outer], [line, line], ls='dashed',
                 color='red', zorder=2, alpha=0.6)

    plt.title("Time: {time:6.3f} Myr".format(time=t/365.25/1e6+100))
    plt.tight_layout()

    plt.savefig(output_folder +
                "frame_{idx:04d}.jpg".format(idx=idx), dpi=200)
    print("Saved! frame: {idx:04d}".format(idx=idx))
    plt.close()
    idx += 1
