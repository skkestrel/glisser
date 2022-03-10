import numpy as np
import matplotlib.pyplot as plt
import struct


file = open("mercury_info.out", 'r')
lines = file.readlines()
discard_time = []
for line in lines:
    seg = line.split()
    discard_time.append(float(seg[-2]))

discard_time = np.array(discard_time)
discard_time = np.sort(discard_time)
print(discard_time)


np.savetxt("discard_mercury_UNR.out", np.column_stack([discard_time, range(1,discard_time.size+1)]))