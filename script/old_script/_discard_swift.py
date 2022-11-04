import numpy as np
import matplotlib.pyplot as plt
import struct

file = "swift_results.txt"


data = np.loadtxt(file)

discard_time = data[:,1]
discard_time = np.sort(discard_time)
print(discard_time)


np.savetxt("discard_swift_avignon_3.out", np.column_stack([discard_time, range(1,discard_time.size+1)]))