import numpy as np
import matplotlib.pyplot as plt
import os
import struct

folder = "/home/yhuang/GLISSER/glisser/"
target_folder = "fast/cold_output/out-ivp-ifree-low-N1/"
# files = [filename1]
                                    
filename = target_folder + "state.out"
output_folder = folder + target_folder + "reoutput/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

endtime = 1460000000000
midtime = 36525000000

num = 0
with open(folder + filename, 'r') as f:
    planet_count = int(f.readline())
    skip = (planet_count - 1) * 4 + 2
    # print(skip)
    content = f.readlines()[skip+2::3]
    idx = []
    lifetime = []
    for line in content:
        state = int(line.split()[1])
        idx.append(int(line.split()[0]))
        if not state:
            lifetime.append(midtime)
        else:
            lifetime.append(int(line.split()[2]))

num = len(idx)
Z = [x for _,x in sorted(zip(idx, lifetime))]


target_folder = "fast/cold_output/out-ivp-ifree-low/"
filename = target_folder + "state.out"
num = 0
with open(folder + filename, 'r') as f:
    planet_count = int(f.readline())
    skip = (planet_count - 1) * 4 + 2
    # print(skip)
    content = f.readlines()[skip+2::3]
    idx = []
    lifetime = []
    for line in content:
        state = int(line.split()[1])
        idx.append(int(line.split()[0]))
        if not state:
            lifetime.append(endtime)
        else:
            lifetime.append(int(line.split()[2]))

num = len(idx)
Z_ref = [x for _,x in sorted(zip(idx, lifetime))]


stable_count = []
count = 0
a = 0
b = 0
for i, t, t_ref in zip(np.arange(1, num+1), Z, Z_ref):
    count += 1
    if t_ref == endtime:
        a += 1
    if t >= midtime:
        b += 1
    if not (count%10):
        if b == 0:
            out = 0
        else:
            out = np.clip(round(a/b*10),0,10)
        stable_count.append(out)
        print(a, b, out)
        a = 0
        b = 0


print(len(stable_count))

output = open(output_folder + "lifetime.txt", "w")
count = 0
for i, t in zip(np.arange(1, num+1), Z):
    count += 1
    idx = (count-1) // 10
    out = count%10
    if out == 0:
        out = 10
    output.write("{0:d} {1:d} {2:d} {3:d}\n".format(i, t, out, stable_count[idx]))
    
output.close()

        
