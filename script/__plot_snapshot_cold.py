import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import struct
import os
import pandas as pd
import matplotlib
import orbitplotkit as opk
import matplotlib.lines as mlines

font = {'size': 13}
matplotlib.rc('font', **font)
light_blue = '#5ed2f1'
dark_blue = '#0084c6'
yellow = '#ffb545'
red = '#ff2521'


# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Cold/"
output_folder = folder + "pics/snapshots/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# time = [0, 182000000000]
# time = [6275000000]
time = np.arange(0, 36600000001, 100000000)
idx = 0
time = [0, 36600000000]


a0 = 41
a1 = 48
q0 = 28
q1 = 48
fig, [ax1, ax2] = plt.subplots(2, figsize=(4.5, 4.2), sharex=True)
xx = np.linspace(40,50,1000)
yy = 1 - 30/xx
ax1.plot(xx, yy, ls='dashed', c='gray', lw=1, alpha=0.8)
for t, size, color in zip(time, [1.5, 3], ['C1', 'C0']):
    pl_txt = "reoutput/ifree-snapshots/planets_{0:d}.txt".format(t)
    df_pl = opk.snapshot2df(folder+pl_txt)
    # df_pl['Omega'] = np.deg2rad(df_pl['Omega'])
    pl_num = len(df_pl.index)
    pa_txt = "reoutput/ifree-snapshots/particles_{0:d}.txt".format(t)
    df_pa = opk.snapshot2df(folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'ifree'])
    if t == 0:
        df_pa0 = df_pa
    df_pa['a0'] = df_pa0['a']
    df_pa['ifree'] = np.rad2deg(df_pa['ifree'])
    df_pa['color'] = opk.BLUE
    df_pa.loc[df_pa['a0'] < 43.5, 'color'] = opk.RED


    # ax1.scatter(df_pl['a'], df_pl['q'], alpha=1, s=150,
    #             color=df_pl['color'], marker='x', lw=3, zorder=2)
    # ax2.scatter(df_pl['a'], df_pl['ifree'], alpha=1,
    #             s=150, color=df_pl['color'], marker='x', lw=3, zorder=2)


    par_size = size
    if t== 0:
        ax1.scatter(df_pa['a'], df_pa['e'], alpha=size/5,
                    s=par_size, edgecolors='none', facecolors=color, rasterized=True)
        ax2.scatter(df_pa['a'], df_pa['ifree'], alpha=size/5,
                    s=par_size, edgecolors='none', facecolors=color, rasterized=True)
    else:
        df_pa = df_pa
        ax1.scatter(df_pa['a'], df_pa['e'], alpha=size/5,
                    s=par_size, edgecolors='none', facecolors=df_pa['color'], rasterized=True)
        ax2.scatter(df_pa['a'], df_pa['ifree'], alpha=size/5,
                    s=par_size, edgecolors='none', facecolors=df_pa['color'], rasterized=True)

    # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
    # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

    aN = df_pl['a'].iloc[-2]
    x = np.linspace(a0, a1, 1000)
    y = np.linspace(a0, a1*10000, 1000)
    # ax1.plot(x, x, color='grey', linestyle='dashed', lw=1)
    # ax1.fill_between(x, x, y, facecolor='gray', alpha=0.5)

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
    idx += 1


blue_dot = mlines.Line2D([], [], color=opk.BLUE, marker='.', linestyle='None',
                                 markersize=8, label='T = 100 Myr')
red_dot = mlines.Line2D([], [], color=opk.RED, marker='.', linestyle='None',
                                 markersize=8, label='T = 0')
# ax1.legend(loc='lower left', fontsize=10,
#                        handles=[red_dot, blue_dot]).set_zorder(100)
ax2.set_xlabel("a (au)")
ax1.set_ylabel("e")
ax2.set_ylabel(r"$I_{free}$ (deg)")

ax1.set_xlim(a0, a1)
ax1.set_ylim(5e-4, 0.5)
ax1.set_yscale('log')

ax2.set_xlim(a0, a1)
ax2.set_ylim(1e-2, 50)
ax2.set_yscale("log")

a_N = aN

# ax1.set_title("Time: {time:6.3f} Myr".format(time=t/365.25/1e6))
fig.tight_layout()
fig.subplots_adjust(hspace=0)
plt.savefig(output_folder+"frame_free_{idx:04d}.jpg".format(idx=idx), dpi=300)
plt.savefig(output_folder+"frame_free_{idx:04d}.pdf".format(idx=idx), dpi=300)
print("Saved! frame: {idx:04d}".format(idx=idx))
plt.close()

