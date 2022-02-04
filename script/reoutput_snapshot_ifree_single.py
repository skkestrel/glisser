import numpy as np
import matplotlib.pyplot as plt
import os
import struct
import rebound
import orbitplotkit as opk
import invariant_plane as ivp

folder = "/home/yhuang/GLISSER/glisser/"
target_folder = "fast/cold_output/out-ivp-ifree-low/"
# files = [filename1]

filename = target_folder + "tracks/track.0.out"
output_folder = folder + target_folder + "reoutput/ifree-snapshots/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

num = 0
# interval = 1
start = 0
end = 1460000000000
output_t = np.arange(start, end, 20000000000)
# print(output_t)
output_t = [40000000000]
with open(folder + filename, 'rb') as f:
    read = f.read(24)
    while len(read) == 24:
        # isOutput = (num%interval == 0)
        time, solar_mass, npl = struct.unpack('=2dQ', read)
        isOutput = (time in output_t)
        # print(time, npl)
        if isOutput:
            output = open(output_folder + "planets_{0:d}.txt".format(int(time)), "w")
        mlist, alist, elist, inclist, Omegalist, omegalist = [], [], [], [], [], []
        for i in range(npl):
            pid, pl_mass, a, e, I, O, o, F = struct.unpack('=I7f', f.read(32))
            mlist.append(pl_mass)
            alist.append(a)
            elist.append(e)
            inclist.append(I)
            Omegalist.append(O)
            omegalist.append(o)
            if isOutput:
                output.write("{id:d} {a:.6f} {e:.6f} {I:.6f} {O:.6f} {o:.6f} {F:.6f}\n"
                            .format(id=pid,a=a,e=e,I=I, O=O, o=o, F=F))
        if isOutput:
            output.close()
            output = open(output_folder + "particles_{0:d}.txt".format(int(time)), "w")
            temp = ivp.InvariantPlane(4, np.array(mlist), np.array(alist), np.array(elist), np.array(inclist), np.array(Omegalist), np.array(omegalist))
        npa, = struct.unpack('=Q', f.read(8))

        for i in range(npa):
            pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
            if isOutput:
                I_free, Omega_free = temp.I_Omega_free_da(a, e, I, O, o)
                print("{0} {1:.4f} {2:.4f} {3:.4f} {4:.4f}".format(pid, a, e, np.rad2deg(I), np.rad2deg(I_free)))
                output.write("{id:d} {a:.6f} {e:.6f} {I:.6f} {O:.6f} {o:.6f} {F:.6f} {If:.6f}\n"
                            .format(id=pid,a=a,e=e,I=I, O=O, o=o, F=F, If=I_free))
        if isOutput:
            output.close()
        read = f.read(24)
        num += 1
