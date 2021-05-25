import numpy as np
import matplotlib.pyplot as plt
import os
import struct

folder = "/home/yhuang/GLISSER/glisser/"
target_folder = "examples/out-JSUNR-1m-80000/"
# files = [filename1]
                                    
filename = target_folder + "tracks/track.0.out"
output_folder = folder + target_folder + "reoutput/snapshots/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.rmdir(output_folder)
    os.makedirs(output_folder)

num = 0
interval = 1
with open(folder + filename, 'rb') as f:
    read = f.read(24)
    while len(read) == 24:
        isOutput = (num%interval == 0)
        time, solar_mass, npl = struct.unpack('=2dQ', read)
        # print(time, npl)
        a1, e1, I1, O1, o1= [], [], [], [], []
        if isOutput:
            output = open(output_folder + "planets_{0:d}.txt".format(int(time)), "w")
        for i in range(npl):
            pid, pl_mass, a, e, I, O, o, F = struct.unpack('=I7f', f.read(32))
            if isOutput:
                output.write("{id:d} {a:.4f} {e:.4f} {I:.4f}\n".format(id=pid,a=a,e=e,I=I))
        if isOutput:
            output.close()
            output = open(output_folder + "particles_{0:d}.txt".format(int(time)), "w")
        npa, = struct.unpack('=Q', f.read(8))
        for i in range(npa):
            pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
            if isOutput:
                output.write("{id:d} {a:.4f} {e:.4f} {I:.4f}\n".format(id=pid,a=a,e=e,I=I))
        if isOutput:
            output.close()
        read = f.read(24)
        num += 1
