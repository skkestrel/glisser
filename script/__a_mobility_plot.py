import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
import struct
import os
import sys
import pandas as pd
import matplotlib
import orbitplotkit as opk
import random

font = {'size': 18}
matplotlib.rc('font', **font)


folder1 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Synthetic-4Gyr-filter/"
output_folder = folder1
label = "GLISSER"


def plot(time):
    qcut = 30
    inner = 50
    outer = 250
    width = 2

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    for t in time:
        pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
        df_pl = opk.snapshot2df(folder1+pl_txt)
        # print(df_pl)
        aN = df_pl.iloc[3]['a']
        print(aN)

        pa_txt = "_rogue_sim_export_4Gyr_q35_detached.csv"
        df = pd.read_csv(pa_txt).set_index('id')
        df['q'] = df['a'] * (1 - df['e'])
        df['ascnode'] = (df['a'] * (1 - df['e']**2)) / \
            (1+df['e']*np.cos(-df['omega']))
        df['inc'] = np.rad2deg(df['inc'])
        df['inc0'] = np.rad2deg(df['inc0'])
        df['Omega'] = np.rad2deg(df['Omega'])
        df['omega'] = np.rad2deg(df['omega'])

        df['da'] = (df['a'] - df['a0'])
        df = df[df['q'] > qcut]

        # axes[0][0].scatter(df['a'], df['da'])
        # axes[0][0].set(xlabel='a (au)', ylabel='da = a - a0 (au)',
        #         xlim=(49, 101), ylim=(-50, 50))


        axes[0][0].scatter(df['a0'], df['q0'])
        axes[0][1].scatter(df['a'], df['q'])

        axes[0][0].set(xlim=(inner, outer))
        axes[0][1].set(xlim=(inner, outer))



        nbins = 50
        xx = np.linspace(inner, outer)
        yy = opk.aProbFuncNorm(xx, [inner, outer], -1.5)
        axes[1][0].hist(df['a0'], bins=nbins, density=True, histtype='step', range=(inner, outer), log=True)
        axes[1][0].plot(xx, yy, ls='dashed', c=opk.DARK_RED, label='a^-1.5')
        axes[1][0].set(yscale='log')
        axes[1][0].legend()

        xx = np.linspace(inner, outer)
        yy = opk.aProbFuncNorm(xx, [inner, outer], -1.5)
        axes[1][1].hist(df['a'], bins=nbins, density=True, histtype='step', range=(inner, outer), log=True)
        axes[1][1].plot(xx, yy, ls='dashed', c=opk.DARK_RED, label='a^-1.5')
        axes[1][1].set(yscale='log')
        axes[1][1].legend()

        # axes[0][1].hist(delta_P, bins=50, density=True, cumulative=True, color=opk.DARK_RED, histtype='step', alpha=1, label='q > {q} au, q0 < 35 au, detached'.format(q=qcut))
        # axes[0][1].grid(ls='dashed', alpha=0.5, c='gray')
        # axes[0][1].set(xlabel='\Delta P', ylabel='Cumulative fraction')

        fig.tight_layout()
        plt.savefig("a_mobility_qcut={q}.jpg".format(q=qcut), dpi=200)


time = [0]
# plot(time, False, False, False)
# plot(time, False, False, True)
plot(time)
