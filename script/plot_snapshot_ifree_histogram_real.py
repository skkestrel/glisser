import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
import os
import pandas as pd
import matplotlib
import orbitplotkit as opk

# pd.set_option("display.max_rows", 99999)
font = {'weight' : 'bold',
        'size'   : 14}
matplotlib.rc('font', **font)
light_blue = '#5ed2f1'
dark_blue = '#0084c6'
yellow = '#ffb545'
red = '#ff2521'



hist_range = (0, 20)

fig, ax1 = plt.subplots(figsize=(8,6))
n_bins = 100
n, bins, patches = ax1.hist(df_group1['q'], n_bins, density=True, histtype='step',
                           cumulative=True, label='Numerical',range=hist_range)



ax1.set_ylabel("Relative ratio")
ax1.set_xlabel("I (au)")
ax1.grid(True)



df = pd.read_csv("if_final_nan.csv")
df['q'] = df['a']*(1-df['e'])
df_non = df[df['RESO'] == 0].copy()
df_all = df[df['RESO'] >= 0].copy()
df_res = df[df['RESO'] > 0].copy()

cut0 = 0
cut1 = 4.2
cut2 = 9
cut3 = 20
df_non = df[df['RESO'] >= 0].copy()
df_low = df_non[df_non['If_da'].between(cut0, cut1)]
df_real_group1 = df_low[df_low['a'].between(left,middle)]
n, bins, patches = ax1.hist(df_real_group1['q'], n_bins, density=True, histtype='step',
                           cumulative=True, label='Real all',range=hist_range)

df_non = df[df['RESO'] == 0].copy()
df_low = df_non[df_non['If_da'].between(cut0, cut1)]
df_real_group1 = df_low[df_low['a'].between(left,middle)]
n, bins, patches = ax1.hist(df_real_group1['q'], n_bins, density=True, histtype='step',
                           cumulative=True, label='Real classical',range=hist_range)
ax1.legend(loc='center left')
ax1.set_title("Group 1 ({0:.1f} < a < {1:.1f} au)".format(left,middle))

plt.savefig(output_folder+"ifree_cold_low_histogram (group 1)_all.jpg".format(idx=idx),dpi=400)
print("Saved! frame: {idx:04d}".format(idx=idx))
plt.tight_layout()
plt.close()
idx += 1