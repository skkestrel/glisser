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


# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNR-Resonant-E/"
output_folder = folder + "pics/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

df = pd.read_csv(folder + "reoutput/lifetime.txt", delimiter=' ')
df['time'] = df['time']/365.25/1e6
print(df)

fig, ax1 = plt.subplots(figsize=(9, 6))

ax1.plot(df['time'], df['count_all'], label='All')
ax1.plot(df['time'], df['count_det'], label='Datached (q > 38 au)')
ax1.plot(df['time'], df['count_sca'], label='Scattering (q < 38 au)')


ax1.set(xlabel='Time (Myr)', ylabel='Count', xlim=[0,120], ylim=[0,1.1e5])

ax1.grid()
plt.legend()
plt.tight_layout()

plt.savefig("stats.jpg", dpi=200)
print("Saved!")
plt.close()
