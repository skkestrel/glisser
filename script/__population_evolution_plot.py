import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import FormatStrFormatter
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
folder2 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUN-Giant-Synthetic-Ref/"
folder3 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Synthetic-1Gyr/"
folder4 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Synthetic-1Gyr-filter/"
output_folder = folder1 + "/pics"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
aN_bary = 30.068010875914965
qcut = 38


def plot(time1, time2):

    num_tot, num_res, num_det, num_tot_r, num_tot3, num_resN52, num_resN31, num_resN41, num_resN83, num_resN73 = [
    ], [], [], [], [], [], [], [], [], []
    fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(
        8, 5), gridspec_kw={'width_ratios': [2.5, 1]}, sharey=True)
    flag = True
    for t in time1:

        pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
        df_pl = pd.read_csv(
            folder1+pl_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')

        pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
        df_pa_r = pd.read_csv(
            folder2+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')
        df_pa_r['q'] = df_pa_r['a'] * (1 - df_pa_r['e'])

        pa_txt = "reoutput/snapshots_RESO/temp_{0:d}.csv".format(t)
        df_pa = pd.read_csv(folder1+pa_txt, delimiter=',').set_index('id')
        df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])

        # df_highq = df_pa[df_pa['q'] > qcut].copy()
        df_arange = df_pa[df_pa['a'].between(50, 100)].copy()
        df_arange_r = df_pa_r[df_pa_r['a'].between(50, 100)].copy()

        df_arange = df_arange[df_arange['q'] > qcut]
        df_arange_r = df_arange_r[df_arange_r['q'] > qcut]
        # df_highq = df_highq[df_highq['a'].between(50, 100)]
        # print(df_highq)

        num_tot.append(df_arange.shape[0])
        num_tot_r.append(df_arange_r.shape[0])
        num_res.append(df_arange[df_arange['RESO'] > 1].shape[0])

        # aN = df_pl['a'].iloc[-2]
        # klist, jlist, a_rlist = opk.genOuterRes(
        #     aN, 49.5, 101, high1=5, high2=50, order_lim=50)

        num_resN52.append(df_arange[df_arange['RESO'] == 2005].shape[0])
        num_resN31.append(df_arange[df_arange['RESO'] == 1003].shape[0])
        num_resN41.append(df_arange[df_arange['RESO'] == 1004].shape[0])
        num_resN83.append(df_arange[df_arange['RESO'] == 3008].shape[0])
        num_resN73.append(df_arange[df_arange['RESO'] == 3007].shape[0])
        if flag and df_arange[df_arange['RESO'] > 1].shape[0] != 0:
            ts = t
            flag = False
        num_det.append(df_arange[df_arange['RESO'] == 0].shape[0])

    for t in time2:
        pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
        df_pa3 = pd.read_csv(
            folder3+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')
        df_pa3['q'] = df_pa3['a'] * (1 - df_pa3['e'])
        df_arange3 = df_pa3[df_pa3['a'].between(50, 100)].copy()
        df_arange3 = df_arange3[df_arange3['q'] > qcut]
        num_tot3.append(df_arange3.shape[0])

    pa_txt = "reoutput/snapshots/temp_0.csv".format(t)
    df_pa_f = pd.read_csv(folder4+pa_txt, delimiter=',').set_index('id')
    df_pa_f['q'] = df_pa_f['a'] * (1 - df_pa_f['e'])
    df_arangef = df_pa_f[df_pa_f['a'].between(50, 100)].copy()
    df_arangef = df_arangef[df_arangef['q'] > qcut]
    num_tot_f = df_arangef.shape[0]
    num_res_f = df_arangef[df_arangef['RESO'] > 1].shape[0]
    num_det_f = df_arangef[df_arangef['RESO'] == 0].shape[0]

    print(num_tot)
    print(num_tot3)
    print(num_tot_f)
    time1 = time1 / (365.25*1e6)
    time2 = time2 / (365.25*1e6) + 40180000000 / (365.25*1e6)
    ts = ts / (365.25*1e6)
    # print(num_tot, num_res, num_det)

    idx = np.argmax(num_res != 0)
    print(idx)
    total = 1e5
    # ax1.plot(time1[idx:], np.array(num_res[idx:])/total,
    #          label='Resonant', color=opk.ORANGE, lw=1.5)

    # ax2.plot(time, np.array(num_resN52)/np.array(num_res),
    #          label='5:2', color=opk.DARK_RED, lw=1.5, ls='solid')
    # ax2.plot(time, np.array(num_resN31)/np.array(num_res),
    #          label='3:1', color=opk.DARK_RED, lw=1.5, ls='dashed')
    # ax2.plot(time, np.array(num_resN83)/np.array(num_res),
    #          label='8:3', color=opk.DARK_RED, lw=1.5, ls='-.')
    # ax2.plot(time, np.array(num_resN73)/np.array(num_res),
    #          label='7:3', color=opk.DARK_RED, lw=1.5, ls='dotted')
    # ax2.plot(time, np.array(num_resN73)/np.array(num_res),
    #          label='Resonant 7:3', color=opk.DARK_RED, lw=1.5, ls='^')
    # ax1.plot(time1[idx:], np.array(num_det[idx:])/total,
    #          label='Detached', color=opk.BLUE, lw=1.5)
    # ax1.plot([ts, ts], [1, 1e5], color=opk.DARK_RED,
    #          lw=2, zorder=6, ls='dashed')

    ax1.plot(time1, np.array(num_tot)/total, label='Rogue sim: q>38 au',
             color=opk.BLACK, lw=2, zorder=10)
    ax1.plot(time1, np.array(num_tot_r)/total, label='Ref sim: q>38 au',
             color=opk.BLACK, lw=1.5, zorder=10, ls='dashed')

    print(time1)
    time1_end = time1[-1]
    num_tot_end = num_tot[-1]

    time2 = np.insert(time2, 0, time1_end)
    num_tot3 = np.insert(num_tot3, 0, num_tot_end)
    print(time2, num_tot3)
    ax2.plot(time2, np.array(num_tot3)/total,
             color=opk.BLACK, lw=2, zorder=10)
    ax2.scatter(time2[-1], num_tot_f/total, marker='x',
                color=opk.BLACK, lw=2, zorder=10, s=70)
    # ax2.plot([time1_end, time2[-1]], [num_res[-1]/total, num_res_f/total], color=opk.ORANGE, lw=1.5, zorder=10,ls='dotted')
    # ax2.scatter(time2[-1], num_res_f/total, marker='x',
    #             color=opk.ORANGE, lw=2, zorder=10, s=70)

    # ax2.plot([time1_end, time2[-1]], [num_det[-1]/total, num_det_f/total], color=opk.BLUE, lw=1.5, zorder=10,ls='dotted')
    # ax2.scatter(time2[-1], num_det_f/total, marker='x',
    #             color=opk.BLUE, lw=2, zorder=10, s=70)
    ax1.set_xlabel("Time (Myr)")

    # ax2.set_xlabel("Time (Myr)")
    ax1.set_ylabel("Population fractions")

    ax1.set_xlim(0, time1_end)
    ax2.set_xlim(time1_end, 1040)
    per1 = 0.125
    ax1.set_ylim(0, per1)
    yvals = ax1.get_yticks()

    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax1.set_yticks([1e-3, 1e-2, 1e-1])
    # ax1.set_yticklabels(['0.1%', '1%', '10%'])

    ax1.legend(loc='upper left', fontsize=11)
    ax1.annotate('Rogue removal', fontsize=16,
                 xy=(time1[-1], per1), xycoords='data',
                 xytext=(time1[-1]-17.4, per1+0.012), textcoords='data',
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
    # ax2.legend(loc='lower right', fontsize=12, ncol=2)
    # ax1.set_yscale('log')
    # ax2.set_yscale('log')
    ax1.set_yticks(np.array([0, 2, 4, 6, 8, 10, 12])/100)
    ax1.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(opk.percent_fmt))

    # ax2.set_ylim(0, 0.12)
    # ax2.set_ylabel("Resonant population fractions")
    # ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    ax1.grid(c=opk.GRAY, alpha=0.5)
    ax2.grid(c=opk.GRAY, alpha=0.5)

    fig.tight_layout()
    fig.subplots_adjust(wspace=0)
    ax1.xaxis.set_label_coords(0.7, -0.12)
    # ax1.set_xlim(0, 100)

    name = "_population_evolution_only_total.jpg"
    plt.savefig(name, dpi=250)

    print("Saved! frame: {name}".format(name=name))


# print(np.arange(0, 24090000001, 100000000))
plot(np.arange(0, 38340000001, 500000000),
     np.arange(0, 328000000000, 5000000000))
