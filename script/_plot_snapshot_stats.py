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


jpl_data = pd.read_csv("sbdb_query_results_1.19.csv")
jpl_data = jpl_data[jpl_data['q'] > 38]
jpl = jpl_data[jpl_data['a'] > 50].copy()
# print(jpl)

# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Full-Scattering/"
output_folder = folder + "pics/stats/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# time = [0, 182000000000]
# time = [6275000000]
time = np.arange(0, 73100000001, 1000000000)
# print(time.shape[0])
idx = 1
# time = [7300000000]
# time = [0]
par_size = 15

pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(0)
df_pa0 = pd.read_csv(
    folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ')
df_pa0['q'] = df_pa0['a'] * (1 - df_pa0['e'])
df_pa0['inc'] = np.rad2deg(df_pa0['inc'])

nbins = 100
for t in time:
    pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
    df_pl = pd.read_csv(
        folder+pl_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ')
    df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
    df_pl['inc'] = np.rad2deg(df_pl['inc'])
    df_pl['Omega'] = np.rad2deg(df_pl['Omega'])
    df_pl['omega'] = np.rad2deg(df_pl['omega'])
    df_pl['f'] = np.rad2deg(df_pl['f'])

    pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
    df_pa = pd.read_csv(
        folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ')
    df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])
    df_pa['q0'] = df_pa0['q']
    df_pa['inc0'] = df_pa0['inc']
    df_pa['inc'] = np.rad2deg(df_pa['inc'])
    df_pa['Omega'] = np.rad2deg(df_pa['Omega'])
    df_pa['omega'] = np.rad2deg(df_pa['omega'])
    df_pa['f'] = np.rad2deg(df_pa['f'])

    df_pa = df_pa[df_pa['a'].between(50, 200)].copy()

    co_list = [-2.5, -3.5]
    x1, x2 = 50, 200
    xx = np.linspace(50, 200, 1000)
    fig, ax1 = plt.subplots(1, figsize=(9, 6))
    for co, color in zip(co_list, [opk.GREEN, opk.RED]):
        scale = -1/(1+co) * x1**(1+co) + 1/(1+co) * x2**(1+co)

        yy = xx**(co)/scale

        ax1.hist(df_pa['a'], nbins, density=True,
                 color=opk.BLUE, histtype='step')
        ax1.plot(xx, yy, color=color, ls='dashed',
                 label="co = {0:.1f}".format(co))


    # a_N = 30
    # klist, jlist, a_rlist = opk.genOuterRes(
    #     a_N, 50, 125, high1=2, high2=50, order_lim=11)
    # for k, j, a_res in zip(klist, jlist, a_rlist):
    #     label = "{0}/{1}".format(j, k)
    #     ax1.text(a_res-0.3, 0.065, label, fontsize=10, rotation=70)
    #     # print(a)
    #     ax1.plot([a_res, a_res], [0, 300], ls='dashed',
    #              color='gray', zorder=1, alpha=0.6)

    ax1.set_xlabel("a (au)")
    ax1.set_ylabel("Count")
    ax1.set_xlim(45, 205)
    ax1.set_ylim(0, 0.06)
    ax1.legend()
    plt.title("Time: {time:6.3f} Myr".format(time=t/365.25/1e6))
    plt.tight_layout()

    plt.savefig(output_folder +
                "a_stats_{idx:04d}.jpg".format(idx=idx), dpi=200)
    print("Saved! frame: {idx:04d}".format(idx=idx))
    plt.close()
    idx += 1
