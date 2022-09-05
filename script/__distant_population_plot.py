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
output_folder = folder1 + "/pics"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
aN_bary = 30.068010875914965


def plot(time):
    # df_temp = pd.read_csv(
    #     folder1+"classes.txt", names=['id', 'RESO', 'k1', 'k2'], delimiter=' ').set_index('id')

    # time = [0, 21520000000]
    idx = 1
    # time = [0]

    par_size = 3
    qcut = 38
    q_low_bound = 32.5
    q_up_bound = 64
    i0, i1 = 0, 62
    per1 = 0.12
    outer = 101
    inner = 49.5

    lower = 1e-1
    upper = 1.5e1

    nbin = 200
    for t in time:
        fig, ax1 = plt.subplots(1, figsize=(8, 6))

        pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
        df_pl = pd.read_csv(
            folder1+pl_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')

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
        df_pa.loc[df_pa['q'] < qcut, 'alpha'] = 0.2
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
        df_arange = df_arange[df_arange['q'] > qcut]
        # df_highq = df_highq[df_highq['a'].between(50, 100)]
        # print(df_highq)

        ax1.set_ylabel("a (au)")
        ax1.set_ylabel("Population ratio relative to 8:3")

        ax1.set_xlim(inner, outer)
        ax1.set_ylim(lower, upper)
        ax1.set_yscale('log')

        aN = df_pl['a'].iloc[-2]
        a_N = aN
        klist, jlist, a_rlist = opk.genOuterRes(
            a_N, 49.5, 101, high1=5, high2=50, order_lim=24)

        num83 = df_arange[df_arange['RESO'] == 3008].shape[0]
        if num83 == 0:
            continue
        print(num83)
        for k, j, a_res in zip(klist, jlist, a_rlist):
            label = "{0:0d}/{1:0d}".format(j, k)
            ax1.text(a_res, upper*1.2, label,
                     fontsize=k**(-0.4)*12, rotation=60, ha='center')
            RESO = j + 1000*k
            num = df_arange[df_arange['RESO'] == RESO].shape[0]
            ax1.scatter(a_res, num/num83, marker='x',
                        facecolors=opk.DARK_RED, edgecolors='none', s=100)

            ax1.plot([a_res, a_res], [lower, upper], ls='dashed',
                     color=opk.GRAY, zorder=-1, lw=1)

            # print(RESO)

        plt.savefig(output_folder +
                    "/distant_population_{idx:04d}_q>38.jpg".format(idx=idx), dpi=250)
        print("Saved! frame: {idx:04d}".format(idx=idx))
        plt.close()
        idx += 1


print(np.arange(0, 24090000001, 1000000000))
plot(np.arange(0, 24090000001, 1000000000))
