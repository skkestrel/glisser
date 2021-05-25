import numpy as np
import matplotlib.pyplot as plt
import struct

folder = "/home/yhuang/GLISSER/glisser/"
# filename1 = "sockeye_results/rogue_JSUNR_80000_1Gyr.out"
filename2 = "benchmark/rogue-particles-JSUNR-1000/tracks/track.0.out"
output_folder = "/home/yhuang/GLISSER/glisser/Rogue_JSUNR_1000/"
# files = [filename1]
                                    
filename = filename2
label = "GLISSER"

num = 0
interval = 10
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
            pid, pl_mass, a, e, I, O, o, F = struct.unpack('=I7d', f.read(60))
            if isOutput:
                output.write("{id:d} {a:.4f} {e:.4f} {I:.4f}\n".format(id=pid,a=a,e=e,I=I))
        if isOutput:
            output.close()
            output = open(output_folder + "particles_{0:d}.txt".format(int(time)), "w")
        npa, = struct.unpack('=Q', f.read(8))
        for i in range(npa):
            pid, a, e, I, O, o, F = struct.unpack('=I6d', f.read(52))
            if isOutput:
                output.write("{id:d} {a:.4f} {e:.4f} {I:.4f}\n".format(id=pid,a=a,e=e,I=I))
        if isOutput:
            output.close()
        read = f.read(24)
        num += 1
