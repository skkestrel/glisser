import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import struct
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

time = [40000000000]
# time = [6275000000]
# time = np.arange(0,365000001,500000)
idx = 1


pa_txt = "ifree-snapshots/particles_{0:d}.txt".format(time[0])
df_pa0 = pd.read_csv(folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'If'], delimiter=' ').set_index('id')
pa_txt = "ifree-snapshots/particles_{0:d}.txt".format(time[1])
df_pa1 = pd.read_csv(folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'If'], delimiter=' ').set_index('id')

pa_txt = "lifetime.txt"
df_lt = pd.read_csv(folder+pa_txt, names=['id','lifetime'], delimiter=' ').set_index('id')

df_subtract = (np.abs(df_pa1 - df_pa0)).dropna(axis=0)
df_pa0['delta_a'] = ((df_pa0.index-1) % 10-5) * 0.008
df_target = df_pa0


df_target['If_diff'] = df_subtract['If']
df_target['lifetime'] = df_lt['lifetime']/365.25/1e6

fig, axes = plt.subplots(2, figsize=(14,14))
ax1, ax2 = axes[0], axes[1]

df_target['q'] = df_target['a'] * (1 - df_target['e'])

# cut = np.deg2rad(1)
# df_pa2 = df_target[df_target['If_diff'] < cut].copy()
df_pa2 = df_target.copy() 
# print(df_pa2)
df_pa2['a']  = df_pa2['a'] + df_pa2['delta_a']

par_size = 2
cut1 = np.deg2rad(4.2)
cut2 = np.deg2rad(9)
cut3 = np.deg2rad(20)
df_pa3 = df_pa2.copy()
df_pa_stable = df_pa3[df_pa3['lifetime'] > 3999].copy()
df_pa_unstable = df_pa3[df_pa3['lifetime'] < 3999].copy()
print(df_pa_unstable)


ax1.scatter(df_pa_unstable['a'], df_pa_unstable['q'], c=df_pa_unstable['lifetime'], norm=matplotlib.colors.LogNorm(), vmin=10, alpha=1, s=par_size, edgecolors='none')
ax2.scatter(df_pa_unstable['a'], np.rad2deg(df_pa_unstable['If']), c=df_pa_unstable['lifetime'], norm=matplotlib.colors.LogNorm(), vmin=10, alpha=1, s=par_size, edgecolors='none')

ax1.scatter(df_pa_stable['a'], df_pa_stable['q'], alpha=1, s=par_size, edgecolors='none', facecolors='black')
ax2.scatter(df_pa_stable['a'], np.rad2deg(df_pa_stable['If']), alpha=1, s=par_size, edgecolors='none', facecolors='black')

fig.colorbar(plt.gca().get_children()[0], ax=axes.ravel().tolist())

# ax1.scatter(df_pa2['a'][df_pa2['If'] < cut1], df_pa2['q'][df_pa2['If'] < cut1], alpha=1, s=par_size, edgecolors='none', facecolors='C0')
# ax2.scatter(df_pa2['a'][df_pa2['If'] < cut1], np.rad2deg(df_pa2['If'][df_pa2['If'] < cut1]), alpha=1, s=par_size, edgecolors='none', facecolors='C0')

# ax1.scatter(df_pa2['a'][df_pa2['If'].between(cut1, cut2)], df_pa2['q'][df_pa2['If'].between(cut1, cut2)], alpha=1, s=par_size, edgecolors='none', facecolors='C0')
# ax2.scatter(df_pa2['a'][df_pa2['If'].between(cut1, cut2)], np.rad2deg(df_pa2['If'][df_pa2['If'].between(cut1, cut2)]), alpha=1, s=par_size, edgecolors='none', facecolors='C0')

# ax1.scatter(df_pa2['a'][df_pa2['If'].between(cut2, cut3)], df_pa2['q'][df_pa2['If'].between(cut2, cut3)], alpha=0.4, s=par_size, edgecolors='none', facecolors='C0')
# ax2.scatter(df_pa2['a'][df_pa2['If'].between(cut2, cut3)], np.rad2deg(df_pa2['If'][df_pa2['If'].between(cut2, cut3)]), alpha=0.4, s=par_size, edgecolors='none', facecolors='C0')

# ax1.scatter(df_pa2['a'][df_pa2['If'] > cut3], df_pa2['q'][df_pa2['If'] > cut3], alpha=0.4, s=par_size, edgecolors='none', facecolors='C4')
# ax2.scatter(df_pa2['a'][df_pa2['If'] > cut3], np.rad2deg(df_pa2['If'][df_pa2['If'] > cut3]), alpha=0.4, s=par_size, edgecolors='none', facecolors='C4')

# ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
# ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

aN =  30
outer = 50
inner = 30
x = np.linspace(inner,outer,1000)
y = np.linspace(inner,outer*10000,1000)
ax1.plot(x, x, color='grey', linestyle='dashed',lw=1)
ax1.fill_between(x, x, y, facecolor='gray', alpha=0.5)

ax1.plot([30, outer], [38, 38], color=red, alpha=0.5, linestyle='dashed',lw=2)
ax2.plot([30, outer], [4.2, 4.2], color=red, linestyle='dashed',lw=2)
ax2.plot([30, outer], [9, 9], color=red, linestyle='dashed',lw=2)

# ax1.plot([aN, outer], [aN, aN], color=red, linestyle='dashed',lw=2)
# ax1.plot([38, outer], [38, 38], color=red, alpha=0.5, linestyle='dashed',lw=2)
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
ax2.set_ylabel("If (deg)")

# ax1.set_title("UNSTABLE\n\n", color='red')
# ax1.set_yscale("log")
# ax2.set_xscale("log")
top_boarder = 42
ax1.set_xlim(39.5,48.5)
ax1.set_ylim(36,top_boarder)

ax2.set_xlim(39.5,49)
ax2.set_ylim(0,5.5)
# ax2.set_yscale("log")

a_N = aN
a_32 = opk.resonanceLocationBYA(a_N, 2, 3)
a_85 = opk.resonanceLocationBYA(a_N, 5, 8)
a_53 = opk.resonanceLocationBYA(a_N, 3, 5)
a_127 = opk.resonanceLocationBYA(a_N, 7, 12)
a_74 = opk.resonanceLocationBYA(a_N, 4, 7)
a_95 = opk.resonanceLocationBYA(a_N, 5, 9)
a_21 = opk.resonanceLocationBYA(a_N, 1, 2)

for a, label in zip([a_32, a_85, a_53, a_127, a_74, a_95, a_21],
                [" 3/2", " 8/5", " 5/3", "12/7", " 7/4", " 9/5", " 2/1"]):
        ax1.annotate(label, fontsize=14,
                        xy=(a, top_boarder), xycoords='data',
                        xytext=(a-0.1, top_boarder+0.8), textcoords='data',
                        arrowprops=dict(arrowstyle="-",
                                        connectionstyle="arc3"),
                        )



# data = np.loadtxt("_nu8_inc=2_new.txt", usecols=[0,1,2])
# a, q = data[:, 0], data[:, 2]
# ax1.plot(a, q, color=red, linestyle='dashed', alpha=0.8, lw=1.6, zorder=-1)

# data = np.loadtxt("_nu8_inc=8_new.txt", usecols=[0,1,2])
# a, q = data[:, 0], data[:, 2]
# ax1.plot(a, q, color=red, linestyle='dashed', alpha=0.8, lw=1.6, zorder=-1)

df = pd.read_csv("if_final_nan.csv")
df['q'] = df['a']*(1-df['e'])
df_non = df[df['RESO'] == 0].copy()
df_res = df[df['RESO'] > 0].copy()

cut0 = 0
cut1 = 4.2
cut2 = 9
cut3 = 20

s = 15
ax1.scatter(df_non['a'][df_non['If_da'].between(cut0, cut1)], df_non['q'][df_non['If_da'].between(cut0, cut1)], facecolors=red, marker='x',
        edgecolors='none', label=r'Classical $I_{free} = (0, 4.2)^\circ$', s=s, alpha=0.6)
# ax2.scatter(df_non['a'][df_non['If_da'].between(cut0, cut1)], df_non['If_da'][df_non['If_da'].between(cut0, cut1)], facecolors=red, marker='x',
#         edgecolors='none', label=r'Classical $I_{free} = (0, 4.2)^\circ$', s=s, alpha=0.6)
# ax1.scatter(df_non['a'][df_non['If_da'].between(cut1, cut2)], df_non['q'][df_non['If_da'].between(cut1, cut2)], facecolors='purple', marker='x',
#         edgecolors='none', label=r'Classical $I_{free} = (4.2, 9)^\circ$', s=s, alpha=0.8)
# ax2.scatter(df_non['a'][df_non['If_da'].between(cut1, cut2)], df_non['If_da'][df_non['If_da'].between(cut1, cut2)], facecolors='purple', marker='x',
#         edgecolors='none', label=r'Classical $I_{free} = (4.2, 9)^\circ$', s=s, alpha=0.8)
# ax1.scatter(df_non['a'][df_non['If_da'].between(cut2, cut3)], df_non['q'][df_non['If_da'].between(cut2, cut3)], facecolors='magenta', marker='x',
#         edgecolors='none', label=r'Classical $I_{free} = (9, 20)^\circ$', s=s, alpha=0.8)
# ax2.scatter(df_non['a'][df_non['If_da'].between(cut2, cut3)], df_non['If_da'][df_non['If_da'].between(cut2, cut3)], facecolors='magenta', marker='x',
#         edgecolors='none', label=r'Classical $I_{free} = (9, 20)^\circ$', s=s, alpha=0.8)


# cut = 4
# ax1.scatter(df_non['a'][df_non['If_da'] < cut], df_non['q'][df_non['If_da'] < cut], facecolors=red,
#         edgecolors='none', label=r'Classical $I_{free} < 4^\circ$', s=s, alpha=0.8)
# ax2.scatter(df_non['a'][df_non['If_da'] < cut], df_non['If_da'][df_non['If_da'] < cut], facecolors=red,
#         edgecolors='none', label=r'Classical $I_{free} < 4^\circ$', s=s, alpha=0.8)

ax1.legend(fancybox=True, scatterpoints=3, fontsize=13)


plt.title("Time: {time:6.3f} Myr".format(time=time[1]/365/1e6))

plt.savefig(output_folder+"ifree_cold_overlay_low.jpg".format(idx=idx),dpi=300)
print("Saved! frame: {idx:04d}".format(idx=idx))
plt.tight_layout()
plt.close()
idx += 1