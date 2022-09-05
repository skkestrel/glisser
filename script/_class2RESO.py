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
pd.set_option('display.max_rows', 30)

# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Giant-Synthetic/"
folder2 = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Synthetic-4Gyr-filter/"
snapshots_RESO_folder = folder + "reoutput/snapshots_RESO/"
snapshots_folder = folder2 + "reoutput/snapshots/"
label = "GLISSER"


df0 = pd.read_csv(snapshots_RESO_folder + "temp_0.csv",
                  delimiter=',').set_index('id')

df = pd.read_csv(snapshots_folder + "particles_0.txt",
                 names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ').set_index('id')

df_class = pd.read_csv(folder2+"classes.txt",
                       names=['id', 'RESO', 'j', 'k'], delimiter=' ').set_index('id').astype(int)

df_class['new_RESO'] = df_class['j']*1000+df_class['k']
df_class.loc[df_class['RESO'] > 0, 'RESO'] = df_class['new_RESO']

df['RESO'] = df_class['RESO']
df['RESO'] = df['RESO'].astype('Int64')
df['q'] = df['a']*(1-df['e'])
df['q0'] = df0['a']*(1-df0['e'])
df['a0'] = df0['a']
df['inc0'] = df0['inc']

df = df[df['q0'] < 35]
df = df[df['a'].between(48, 500)]
print(df)
print(len(df[df['RESO'] == 0].index))
print(len(df[df['RESO'] > 0].index))
print(len(df[df['RESO'] < 0].index))
# df = df[df['RESO'] == 0]

# print(df)

# df.to_csv("rogue_sim_export_4Gyr_detached.csv")
