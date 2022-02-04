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


folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
output_folder = folder + "pics/ifree-snapshots/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

time = [0, 40000000000, 1460000000000]
# time = [6275000000]
# time = np.arange(0,365000001,500000)
idx = 1
# qstep = 0.06
astep = 0.02
left = 42.4
middle = 43.7
right = 44.5
rr = 47.1

left = left
middle = middle

pa_txt = "ifree-snapshots/particles_{0:d}.txt".format(time[0])
df_pa0 = pd.read_csv(folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'If'], delimiter=' ').set_index('id')
pa_txt = "ifree-snapshots/particles_{0:d}.txt".format(time[1])
df_pa1 = pd.read_csv(folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'If'], delimiter=' ').set_index('id')
pa_txt = "ifree-snapshots/particles_{0:d}.txt".format(time[2])
df_pa2 = pd.read_csv(folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'If'], delimiter=' ').set_index('id')

pa_txt = "lifetime.txt"
df_lt = pd.read_csv(folder+pa_txt, names=['id','lifetime', 'idx', 'stable_count'], delimiter=' ').set_index('id')
print(df_lt)

df_target = df_pa0
df_target['lifetime'] = df_lt['lifetime']/365.25/1e6
df_target['idx'] = df_lt['idx']
df_target['stable_count'] = df_lt['stable_count']


df_target['q'] = df_target['a'] * (1 - df_target['e'])

# df_pa2 = df_target[df_target['If_diff'] < cut].copy()
df_pa2 = df_target.copy()
# print(df_pa2)

par_size = 2
cut1 = np.deg2rad(4.2)
cut2 = np.deg2rad(9)
cut3 = np.deg2rad(20)
df_pa3 = df_pa2.copy()

fig, ax1 = plt.subplots(figsize=(8,6))
hist_range = (36, middle)
n_bins = 300


df_pa_stable = df_pa3[df_pa3['lifetime'] > 100].copy()
df_group1 = df_pa_stable[df_pa_stable['a'].between(left-astep, middle+astep)]
# print(df_group1)
n, bins, patches = ax1.hist(df_group1['q'], n_bins, density=True, histtype='step',
                           cumulative=True, label='100 Myr',range=hist_range)


df_pa_stable = df_pa3[df_pa3['lifetime'] > 100].copy()
df_group1 = df_pa_stable[df_pa_stable['a'].between(left-astep, middle+astep)]
df_group1 = df_group1[df_group1['stable_count'] > df_group1['idx']]
print(df_group1)
n, bins, patches = ax1.hist(df_group1['q'], n_bins, density=True, histtype='step',
                           cumulative=True, label='synthetic',range=hist_range)

df_pa_stable = df_pa3[df_pa3['lifetime'] > 3997].copy()
df_group1 = df_pa_stable[df_pa_stable['a'].between(left-astep, middle+astep)]
n, bins, patches = ax1.hist(df_group1['q'], n_bins, density=True, histtype='step',
                           cumulative=True, label='4 Gyr',range=hist_range)



ax1.set_ylabel("Cumulative ratio")
ax1.set_xlabel("q (au)")
ax1.grid(True)



# df = pd.read_csv("if_final_nan.csv")
# df['q'] = df['a']*(1-df['e'])
# df_non = df[df['RESO'] == 0].copy()
# df_all = df[df['RESO'] >= 0].copy()
# df_res = df[df['RESO'] > 0].copy()

# cut0 = 0
# cut1 = 4.2
# cut2 = 9
# cut3 = 20
# df_non = df[df['RESO'] >= 0].copy()
# df_low = df_non[df_non['If_da'].between(cut0, cut1)]
# df_real_group1 = df_low[df_low['a'].between(left,middle)]
# n, bins, patches = ax1.hist(df_real_group1['q'], n_bins, density=True, histtype='step',
#                            cumulative=True, label='Real all',range=hist_range)

# df_non = df[df['RESO'] == 0].copy()
# df_low = df_non[df_non['If_da'].between(cut0, cut1)]
# df_real_group1 = df_low[df_low['a'].between(left,middle)]
# n, bins, patches = ax1.hist(df_real_group1['q'], n_bins, density=True, histtype='step',
#                            cumulative=True, label='Real classical',range=hist_range)
ax1.legend(loc='center left')
ax1.set_title("Group 1 ({0:.1f} < a < {1:.1f} au)".format(left,middle))

plt.savefig(output_folder+"ifree_cold_low_histogram (group 3)_compare.jpg".format(idx=idx),dpi=400)
print("Saved! frame: {idx:04d}".format(idx=idx))
plt.tight_layout()
plt.close()
idx += 1