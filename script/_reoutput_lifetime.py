import numpy as np
import matplotlib.pyplot as plt
import os
import struct
import pandas as pd

folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNR-Resonant-E/"
# files = [filename1]

output_folder = folder + "reoutput/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

time = [0]
time = np.arange(0, 36520000001, 20000000)
# idx = 1

output = open(output_folder + "lifetime.txt", "w")
output.write("time count_all count_det count_sca\n")
for t in time:
    pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
    df_pa = pd.read_csv(
        folder+pa_txt, names=['id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ').set_index('id')
    df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])
    count_all = df_pa['a'].count()

    df_det = df_pa['a'][df_pa['q'] > 38].copy()
    df_sca = df_pa['a'][df_pa['q'] < 38].copy()
    count_det = df_det.count()
    count_sca = df_sca.count()
    
    print(t, count_all, count_det, count_sca)
    output.write("{0:d} {1:d} {2:d} {3:d}\n".format(t, count_all, count_det, count_sca))

output.close()

