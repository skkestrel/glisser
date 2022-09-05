import numpy as np
import matplotlib.pyplot as plt
import os
import struct

folder = "/home/yhuang/GLISSER/glisser/"
target_folder = "fast/rogue_output/out-JSUNT-Synthetic-4Gyr-filter/"
# files = [filename1]

filename = target_folder + "tracks/track.0.out"
output_folder = folder + target_folder + "reoutput/hist/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

interval = 1
chunk = 2000
counter = 0
timelist = []
with open(folder + filename, 'rb') as f:
    print("Reading track: {0}".format(folder + filename))
    read = f.read(24)
    while len(read) == 24:
        isOutput = (counter % interval == 0)
        time, solar_mass, npl = struct.unpack('=2dQ', read)
        if isOutput:
            timelist.append(time)
        if counter == 0:
            pid1 = {}
            a1, e1, I1, O1, o1, F1 = np.empty((chunk, npl)), np.empty((chunk, npl)), np.empty((chunk, npl)),  \
                np.empty((chunk, npl)), np.empty(
                    (chunk, npl)), np.empty((chunk, npl))
            a1[:], e1[:], I1[:], O1[:], o1[:], F1[:
                                                  ] = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        for i in range(npl):
            if isOutput:
                inner_counter = counter // interval
                pid, pl_mass, a, e, I, O, o, F \
                    = struct.unpack('=I7f', f.read(32))
                if counter == 0:
                    pid1[pid] = i
                idx = pid1[pid]
                a1[inner_counter, idx] = a
                e1[inner_counter, idx] = e
                I1[inner_counter, idx] = I
                O1[inner_counter, idx] = O
                o1[inner_counter, idx] = o
                F1[inner_counter, idx] = F
            else:
                f.read(32)
            # output = open(output_folder + "pl_{0:d}.txt".format(pid), "a")
            # output.write("{time:d} {a:.6f} {e:.6f} {I:.6f} {O:.6f} {o:.6f} {F:.6f}\n"
            #              .format(time=int(time), a=a, e=e, I=I, O=O, o=o, F=F))
            # output.close()
        npa, = struct.unpack('=Q', f.read(8))
        print(int(time), npl, npa)

        if counter == 0:
            pid2 = {}
            a2, e2, I2, O2, o2, F2 = np.empty((chunk, npa)), np.empty((chunk, npa)), np.empty((chunk, npa)),  \
                np.empty((chunk, npa)), np.empty(
                    (chunk, npa)), np.empty((chunk, npa))
            a2[:], e2[:], I2[:], O2[:], o2[:], F2[:] = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        for i in range(npa):
            if isOutput:
                inner_counter = counter // interval
                pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
                if counter == 0:
                    pid2[pid] = i
                idx = pid2[pid]
                a2[inner_counter, idx] = a
                e2[inner_counter, idx] = e
                I2[inner_counter, idx] = I
                O2[inner_counter, idx] = O
                o2[inner_counter, idx] = o
                F2[inner_counter, idx] = F
            else:
                f.read(28)
            # if isOutput:
                # output = open(output_folder + "pa_{0:d}.txt".format(pid), "a")
                # output.write("{time:d} {a:.6f} {e:.6f} {I:.6f} {O:.6f} {o:.6f} {F:.6f}\n"
                #              .format(time=int(time), a=a, e=e, I=I, O=O, o=o, F=F))
                # output.close()
        read = f.read(24)
        counter += 1

    print("Output to: {0}".format(output_folder))
    for i, idx in enumerate(pid1):
        print("Output: pl_{0:d}.txt".format(idx))
        output = open(output_folder + "pl_{0:d}.txt".format(idx), "w+")
        for time, a, e, I, O, o, F in zip(timelist, a1[:counter, i], e1[:counter, i], I1[:counter, i], O1[:counter, i], o1[:counter, i], F1[:counter, i]):
            if np.isnan(a):
                break
            output.write("{time:d} {a:.6f} {e:.6f} {I:.6f} {O:.6f} {o:.6f} {F:.6f}\n"
                         .format(time=int(time), a=a, e=e, I=I, O=O, o=o, F=F))
        output.close()
    for i, idx in enumerate(pid2):
        print("Output: pa_{0:d}.txt".format(idx))
        output = open(output_folder + "pa_{0:d}.txt".format(idx), "w+")
        for time, a, e, I, O, o, F in zip(timelist, a2[:counter, i], e2[:counter, i], I2[:counter, i], O2[:counter, i], o2[:counter, i], F2[:counter, i]):
            if np.isnan(a):
                break
            output.write("{time:d} {a:.6f} {e:.6f} {I:.6f} {O:.6f} {o:.6f} {F:.6f}\n"
                         .format(time=int(time), a=a, e=e, I=I, O=O, o=o, F=F))
        output.close()
