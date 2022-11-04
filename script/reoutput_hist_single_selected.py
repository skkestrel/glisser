import numpy as np
import matplotlib.pyplot as plt
import os
import struct
import pandas as pd
import orbitplotkit as opk

folder = "/home/yhuang/GLISSER/glisser/"
target_folder = "fast/rogue_output/out-JSUNT1-Migration-1EM/"
# files = [filename1]
snapshot_folder = folder + target_folder + "reoutput/snapshots/"

filename = target_folder + "tracks/track.0.out"
output_folder = folder + target_folder + "reoutput/hist/"
label = "GLISSER"

final_snapshot = "planets_54800000000.txt"
df_pl = pd.read_csv(snapshot_folder+final_snapshot, names=[
    'id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ').set_index('id')

final_snapshot = "particles_54800000000.txt"
df_pa = pd.read_csv(snapshot_folder+final_snapshot, names=[
    'id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ').set_index('id')
df_pa['q'] = df_pa['a']*(1-df_pa['e'])
df_q = df_pa[df_pa['q'] > 50]

a_N = df_pl.iloc[3]['a']

# a_41 = opk.resonanceLocationBYA(a_N, 1, 4)
# print(a_N, a_41)
# a_52 = opk.resonanceLocationBYA(a_N, 2, 5)

# hw = 1
# df_a = df_q[df_q['a'].between(a_41 - hw, a_41 + hw)]

idx = np.array(df_q.index)
print(idx)

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def reoutput():
    with open(folder + filename, 'rb') as f:
        read = f.read(24)
        while len(read) == 24:
            time, solar_mass, npl = struct.unpack('=2dQ', read)
            # print(time, npl)
            a1, e1, I1, O1, o1 = [], [], [], [], []
            for i in range(npl):
                pid, pl_mass, a, e, I, O, o, F = struct.unpack('=I7f', f.read(32))
                output = open(output_folder + "pl_{0:d}.txt".format(pid), "a")
                output.write("{time:d} {a:.6f} {e:.6f} {I:.6f} {O:.6f} {o:.6f} {F:.6f}\n"
                            .format(time=int(time), a=a, e=e, I=I, O=O, o=o, F=F))
                output.close()
            npa, = struct.unpack('=Q', f.read(8))
            for i in range(npa):
                pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
                # isOutput = (pid % interval == 0)
                isOutput = (pid in idx)
                if isOutput:
                    print("Writing particle: {0} at {1}".format(pid, time))
                    output = open(output_folder + "pa_{0:d}.txt".format(pid), "a")
                    output.write("{time:d} {a:.6f} {e:.6f} {I:.6f} {O:.6f} {o:.6f} {F:.6f}\n"
                                .format(time=int(time), a=a, e=e, I=I, O=O, o=o, F=F))
                    output.close()
            read = f.read(24)

# reoutput()