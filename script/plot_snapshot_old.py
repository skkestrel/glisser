import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import pandas as pd
import matplotlib
import orbitplotkit as opk

font = {'weight' : 'bold',
        'size'   : 14}
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
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUN-Resonant/reoutput/"
output_folder = folder + "pics/snapshots/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# time = [0, 182000000000]
# time = [6275000000]
time = np.arange(0,36524000001,40000000)
print(time.shape[0])
idx = 1
# time = [0, 36524000000]
par_size = 10

for t in time:
    pl_txt = "snapshots/planets_{0:d}.txt".format(t)
    df_pl = pd.read_csv(folder+pl_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ')
    # print(df_pl)

    pa_txt = "snapshots/particles_{0:d}.txt".format(t)  
    df_pa = pd.read_csv(folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ')
    print(df_pa)

    fig, [ax1, ax2] = plt.subplots(2, figsize=(16,10))
    # df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
    # ax1.scatter(jpl['a'], jpl['q'], alpha=0.5, s=30, color=red, marker='x',lw=2, zorder=2)
    # ax2.scatter(jpl['a'], jpl['i'], alpha=0.5, s=30, color=red, marker='x',lw=2, zorder=2)
    # ax2.scatter(df_pl['a'], np.rad2deg(df_pl['inc']), alpha=0.8, s=100, color=red, marker='x',lw=2, zorder=2)

    df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])
    
    ax1.scatter(df_pa['a'][df_pa['q'] > 38], df_pa['q'][df_pa['q'] > 38], alpha=0.8, s=par_size, edgecolors='none', facecolors='C0')
    ax1.scatter(df_pa['a'][df_pa['q'] < 38], df_pa['q'][df_pa['q'] < 38], alpha=0.15, s=par_size, edgecolors='none', facecolors='C0')
    ax2.scatter(df_pa['a'][df_pa['q'] > 38], np.rad2deg(df_pa['inc'][df_pa['q'] > 38]), alpha=0.8, s=par_size, edgecolors='none', facecolors='C0')
    ax2.scatter(df_pa['a'][df_pa['q'] < 38], np.rad2deg(df_pa['inc'][df_pa['q'] < 38]), alpha=0.15, s=par_size, edgecolors='none', facecolors='C0')
    # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
    # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

    aN =  df_pl['a'].iloc[-1]
    outer = 127
    inner = 48
    x = np.linspace(inner,outer,1000)
    y = np.linspace(inner,outer*10000,1000)
    ax1.plot(x, x, color='grey', linestyle='dashed',lw=1)
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
    ax1.set_xlim(inner_bound,outer_bound)
    ax2.set_xlim(inner_bound,outer_bound)

    q_low_bound = 28
    q_up_bound = 60
    ax1.set_ylim(q_low_bound,q_up_bound)
    ax2.set_ylim(0,60)
    # ax2.set_yscale("log")

    a_N = aN
    a_52 = opk.resonanceLocationBYA(a_N, 2, 5)
    a_72 = opk.resonanceLocationBYA(a_N, 2, 7)
    a_92 = opk.resonanceLocationBYA(a_N, 2, 9)
    a_112 = opk.resonanceLocationBYA(a_N, 2, 11)
    a_132 = opk.resonanceLocationBYA(a_N, 2, 13)
    a_73 = opk.resonanceLocationBYA(a_N, 3, 7)
    a_83 = opk.resonanceLocationBYA(a_N, 3, 8)
    a_103 = opk.resonanceLocationBYA(a_N, 3, 10)
    a_113 = opk.resonanceLocationBYA(a_N, 3, 11)
    a_133 = opk.resonanceLocationBYA(a_N, 3, 13)
    a_143 = opk.resonanceLocationBYA(a_N, 3, 14)
    a_163 = opk.resonanceLocationBYA(a_N, 3, 16)
    a_173 = opk.resonanceLocationBYA(a_N, 3, 17)
    a_193 = opk.resonanceLocationBYA(a_N, 3, 19)
    a_203 = opk.resonanceLocationBYA(a_N, 3, 20)
    a_31 = opk.resonanceLocationBYA(a_N, 1, 3)
    a_41 = opk.resonanceLocationBYA(a_N, 1, 4)
    a_51 = opk.resonanceLocationBYA(a_N, 1, 5)
    a_61 = opk.resonanceLocationBYA(a_N, 1, 6)
    a_71 = opk.resonanceLocationBYA(a_N, 1, 7)

    for a, label in zip([a_52, a_31, a_72, a_41, a_92, a_51, a_112, a_61, a_132, a_71, a_73, a_83, a_103, a_113, a_133, a_143, a_163, a_173, a_193, a_203],
                        [" 5/2", " 3/1", " 7/2", " 4/1", " 9/2", " 5/1", "11/2", " 6/1", "13/2", " 7/1", " 7/3", " 8/3", "10/3", "11/3", "13/3", "14/3", 
                        "16/3", "17/3", "19/3", "20/3"]):
        ax1.text(a-1, q_up_bound+1, label, fontsize=12, rotation=60)
        # print(a)
        ax1.plot([a, a], [0, 100], ls='dashed',color='gray',zorder=1,alpha=0.6)

    hlines = [30, 38]
    for line in hlines:
        ax1.plot([inner, outer], [line, line], ls='dashed',color='red',zorder=2,alpha=0.6)


    plt.title("Time: {time:6.3f} Myr".format(time=t/365/1e6))

    plt.savefig(output_folder+"frame_{idx:04d}.jpg".format(idx=idx),dpi=200)
    print("Saved! frame: {idx:04d}".format(idx=idx))
    plt.close()
    idx += 1 