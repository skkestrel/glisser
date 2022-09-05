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
    qcut = 38

    for t in time:
        pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
        df_pl = pd.read_csv(
            folder1+pl_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')
        df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
        df_pl['inc'] = np.rad2deg(df_pl['inc'])
        df_pl['Omega'] = np.rad2deg(df_pl['Omega'])
        df_pl['omega'] = np.rad2deg(df_pl['omega'])
        df_pl['M'] = np.rad2deg(df_pl['M'])

        pa_txt = "reoutput/snapshots/temp_{0:d}.csv".format(t)
        df_pa = pd.read_csv(folder1+pa_txt, delimiter=',').set_index('id')
        df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])
        df_pa['inc'] = np.rad2deg(df_pa['inc'])

        print("Time: {0:.1f} Myr".format(t/365.25/1e6))

        df_pa = df_pa[df_pa['q0'] < 37]
        df_pa = df_pa[df_pa['a'].between(47, 250)]
        total_num = len(df_pa)
        print("Total: {0:d}".format(total_num))

        det_num = len(df_pa[df_pa['RESO'] == 0])
        sca_num = len(df_pa[df_pa['RESO'] == -1])
        print("Detached: {0:d}".format(det_num))
        print("Scattering: {0:d}".format(sca_num))
        print("Detached/Scattering Ratio: {0:.2f}".format(det_num/sca_num))
        # print(df_pa)

        df_pa['alpha'] = 1
        # df_pa.loc[df_pa['RESO'] == 1, 'alpha'] = 0.5
        df_pa.loc[df_pa['q'] < qcut, 'alpha'] = 0
        # df_pa.loc[df_pa['q'] > 40, 'alpha'] = 1

        df_pa['parsize'] = 1.5
        df_pa.loc[df_pa['q'] > qcut, 'parsize'] = 5
        # df_pa.loc[df_pa['q'] > 40, 'parsize'] = 4

        df_pa['color'] = opk.GRAY

        df_highq = df_pa[df_pa['q'] > qcut].copy()
        df_arange = df_highq[df_highq['a'].between(50, 100)].copy()
        df_used = df_arange[df_arange['RESO'] == 0]
        # df_highq = df_highq[df_highq['a'].between(50, 100)]
        # print(df_highq)


time = [36520000000]
time = [0]
# plot(time, False, False, False)
# plot(time, False, False, True)
plot(time)
