import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
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
aN = 30.11

# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Giant-Synthetic/"
output_folder = folder + "pics/polar/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


def plot(time, idx, plotReal=False, plotCold=400000):

    # idx = 226
    # time = [0]
    aleft = 22
    a0 = 200
    a1 = 1200
    q0 = 22
    q1 = 80
    qcut = 40
    incmax = 55
    for t in time:
        pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
        df_pl = opk.snapshot2df(folder+pl_txt)
        df_pl = df_pl[df_pl.index > 4]
        # df_pl['varpi'] = np.deg2rad(df_pl['varpi'])
        # df_pl['Omega'] = np.deg2rad(df_pl['Omega'])
        pl_num = len(df_pl.index)
        df_pl['color'] = [opk.ORANGE][:pl_num]
        pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
        df_pa = opk.snapshot2df(folder+pa_txt)
        df_pa = df_pa[df_pa.index <= plotCold]
        df_pa = df_pa[df_pa['q'] >= 25]
        # df_pa['varpi'] = np.deg2rad(df_pa['varpi'])
        # df_pa['Omega'] = np.deg2rad(df_pa['Omega'])



        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(7,6))
        fig.suptitle("Time: {time:6.3f} Myr | a = {a:3.1f} au".format(
            time=t/365.25/1e6, a=df_pl['a'][5]), fontsize=15, y=0.95)
        # fig, [ax3, ax4, ax1, ax2] = plt.subplots(2, 2, figsize=(8.5, 6), gridspec_kw={
        #     'height_ratios': [3, 1], 'width_ratios': [2, 1]}, sharex=True)

        # ax3.scatter(df_pl['a'], df_pl['Omega'], alpha=1,
        #             s=150, color=df_pl['color'], marker='x', lw=3, zorder=2)

        par_size = 0.25
        df_pa['alpha'] = 1
        df_pa['size'] = par_size
        df_pa['color'] = opk.DARK_GRAY
        # df_pa.loc[df_pa.index > plotCold, 'color'] = opk.DARK_RED

        df_highq = (df_pa['a'] > 50) & opk.isDetached(
            df_pa['a'], df_pa['q'], aN=aN, qcut=qcut)

        df_pa.loc[df_highq, 'size'] = par_size*20
        df_pa.loc[df_highq, 'color'] = opk.LIGHT_BLUE
        df_pa.loc[df_highq, 'alpha'] = 1
        df_pa = df_pa.sort_values(by=['q'])

        varpi_5, q5 = df_pl['varpi'][5], df_pl['q'][5]
        print(np.rad2deg(varpi_5), df_pl['Omega'][5], df_pl['omega'][5])
        ax.scatter([varpi_5], [q5], alpha=1, s=100,
                    color=df_pl['color'], marker='x', lw=3, zorder=20, rasterized=True)
        ax.scatter(df_pa['varpi'], df_pa['q'], alpha=df_pa['alpha'],
                       s=df_pa['size'], edgecolors='none', facecolors=df_pa['color'], rasterized=True)
        delta_rad = np.deg2rad(15)
        ax.fill_between(np.linspace(varpi_5-delta_rad,varpi_5+delta_rad, 100),0,q1, color=opk.ORANGE, alpha=0.2,zorder=-10)
        ax.set_rmax(q1)
        # ax.set_rmin(10)
        ax.set_rticks([30, 40, 50, 60, 70])
        ax.set_rlabel_position(-22.5)
        ax.grid(True)
        ax.set_ylabel("q (au)")
        ax.yaxis.labelpad = 35
        ax.set_xlabel(r"$\varpi$ (deg)")


        fig.tight_layout()

        plt.savefig(output_folder +
                    "frame_{idx:04d}.jpg".format(idx=idx), dpi=300)
        print("Saved! frame: {idx:04d}".format(idx=idx))
        plt.close()
        idx += 1


idx = 0
time_step = 50000000
t0 = int(idx*time_step)
time = np.arange(t0, 2*36000000001, time_step)
# time = [t0]
# idx = 31500000000/time_step
plot(time, idx, plotReal=False, plotCold=400000)
