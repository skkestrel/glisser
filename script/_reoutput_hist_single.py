import numpy as np
import matplotlib.pyplot as plt
import os
import struct
import pandas as pd

folder = "/home/yhuang/GLISSER/glisser/"
target_folder = "fast/rogue_output/out-JSUNT1-Migration-1EM/"
# files = [filename1]

filename = target_folder + "tracks/track.0.out"
output_folder = folder + target_folder + "reoutput/hist/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# idlist = pd.read_csv(
#     'JSUNT1-Migration-1EM-outer-particles.csv')['id'].to_list()
# print(idlist)
# minqlist = idlist.copy()
# minqlist[:] = [np.inf for _ in minqlist]


def outputHist(planetOnly=False):
    num = 0
    interval = 1
    with open(folder + filename, 'rb') as f:
        print("Reading track: {0}".format(folder + filename))
        print("Output to: {0}".format(output_folder))
        read = f.read(24)
        while len(read) == 24:
            time, solar_mass, npl = struct.unpack('=2dQ', read)
            # print(time, npl)
            a1, e1, I1, O1, o1 = [], [], [], [], []
            for i in range(npl):
                pid, pl_mass, a, e, I, O, o, F = struct.unpack(
                    '=I7f', f.read(32))
                output = open(output_folder + "pl_{0:d}.txt".format(pid), "a")
                output.write("{time:d} {a:.6f} {e:.6f} {I:.6f} {O:.6f} {o:.6f} {F:.6f}\n"
                             .format(time=int(time), a=a, e=e, I=I, O=O, o=o, F=F))
                output.close()
            npa, = struct.unpack('=Q', f.read(8))
            print(int(time), npl, npa)
            for i in range(npa):
                pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
                isOutput = (pid % interval == 0)
                # isOutput = (pid in idlist)
                if isOutput and (not planetOnly):
                    # idx = idlist.index(pid)
                    # if minqlist[idx] > a*(1-e):
                    #     minqlist[idx] = a*(1-e)
                    output = open(output_folder +
                                  "pa_{0:d}.txt".format(pid), "a")
                    output.write("{time:d} {a:.6f} {e:.6f} {I:.6f} {O:.6f} {o:.6f} {F:.6f}\n"
                                 .format(time=int(time), a=a, e=e, I=I, O=O, o=o, F=F))
                    output.close()
            read = f.read(24)
            num += 1

outputHist(planetOnly=True)
# print(minqlist)
# print(np.mean(minqlist), np.min(minqlist), np.max(minqlist))