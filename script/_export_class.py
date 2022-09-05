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
snapshots_RESO_folder = folder + "reoutput/snapshots_RESO/"
label = "GLISSER"


df0 = pd.read_csv(snapshots_RESO_folder  + "temp_0.csv", delimiter=',').set_index('id')

df = pd.read_csv(snapshots_RESO_folder  + "temp_35020000000.csv", delimiter=',').set_index('id')

df['q'] = df['a']*(1-df['e'])
df['q0'] = df0['a']*(1-df0['e'])
df['a0'] = df0['a']
df['inc0'] = df0['inc']

df.to_csv("distant_export_v1.csv")
