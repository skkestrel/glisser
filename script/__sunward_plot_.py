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
from scipy import stats

font = {'size': 18}
matplotlib.rc('font', **font)


# folder1 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Synthetic-4Gyr-filter/"
# output_folder = folder1
# label = "GLISSER"

# df_des = pd.read_csv("_DETACHED_ALL.csv")

def plot(time):
    qcut = 40
    qup = 60
    inner = 50
    outer = 100
    width = 1
    aN = 30.0687

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 9))
    for t in time:
        
        pa_txt = "_DETACHED_ALL.csv"
        df = pd.read_csv(pa_txt)
        # df['q'] = df['a'] * (1 - df['e'])
        # df['ascnode'] = (df['a'] * (1 - df['e']**2)) / \
        #     (1+df['e']*np.cos(-df['omega']))
        # df['inc'] = np.rad2deg(df['inc'])
        # df['inc0'] = np.rad2deg(df['inc0'])
        # df['Omega'] = np.rad2deg(df['Omega'])
        # df['omega'] = np.rad2deg(df['omega'])

        print(len(df.index))
        klist, jlist, a_rlist = opk.genOuterRes(
            aN, inner, outer, high1=2, high2=7, order_lim=8)
        df['fountain'] = False
        df['a_res'] = 0
        df['color'] = opk.BLUE
        for k, j, a_res in zip(klist, jlist, a_rlist):
            print(k, j, a_res)
            isInside = df['a'].between(a_res - width, a_res + width)
            isRight = df['a'].between(a_res, a_res + width)
            # print(isInside)
            df.loc[isInside, 'fountain'] = True
            df.loc[isInside, 'a_res'] = a_res
            df.loc[isRight, 'color'] = opk.DARK_RED
            label = "{0:0d}/{1:0d}".format(j, k)
            ax1.text(a_res, qup-1.5, label,
                     fontsize=k**(-0.4)*18, rotation=60, ha='center')
            # # print(a)
            ax1.plot([a_res, a_res], [qcut-1, 300], ls='dashed',
                     color=opk.GRAY, zorder=-1, lw=1)
            # ax3.plot([a_res, a_res], [0, per1], ls='dashed',
            #          color=opk.GRAY, zorder=-1, lw=1)
        df_fountain = df[df['fountain'] == True]
        df_fountain = df_fountain[df_fountain['q'] > qcut]
        # print(df_fountain)

        right = len(df_fountain[(df_fountain['color'] == opk.DARK_RED)].index)
        all = len(df_fountain.index)

        P_B = stats.binom_test(right, n=all, p=0.5)

        ax1.scatter(df_fountain['a'], df_fountain['q'],s=10, c=df_fountain['color'])
        ax1.set(xlabel='a (au)', ylabel='q (au)', ylim=(qcut-1, qup), title="left={0}, right={1}, P_B={2:.3f}%".format(all-right, right, P_B*100))

        P = df_fountain['a']**(3/2)
        P_res = df_fountain['a_res']**(3/2)
        P_N = aN**(3/2)

        delta_P = (P - P_res)/P_N
        ax2.hist(delta_P, bins=50, density=True, cumulative=True, color=opk.DARK_RED, histtype='step', alpha=1, label='q > {q} au, q0 < 35 au, detached'.format(q=qcut))
        ax2.grid(ls='dashed', alpha=0.5, c='gray')
        ax2.set(xlabel='Delta P', ylabel='Cumulative fraction', xlim=(-0.15,0.15), title="q > {0} au, width = {1} au".format(qcut, width))

        fig.tight_layout()
        plt.savefig("fountain_qcut={qcut}_width={width}_ALL_DETACHED.jpg".format(qcut=qcut, width=width), dpi=200)


time = [0]
# plot(time, False, False, False)
# plot(time, False, False, True)
plot(time)
