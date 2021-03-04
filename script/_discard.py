import numpy as np
import matplotlib.pyplot as plt
import struct

filename1 = "discard_opio_onlyN_5.out"
filename2 = "discard_mercury_onlyN_5.out"
# filename3 = "discard_glisser_sockeye_100_4.out"
# filename4 = "discard_glisser_sockeye_10_4.out"
# filename5 = "discard_mercury_3.out"
# filename6 = "discard_swift_avignon_3.out"

# filename2 = "benchmark/kuiper-particles-out-1000-read/tracks/track.0.out"

files = [filename1, filename2]
# files = [filename1, filename5]

fig, ax = plt.subplots(figsize=(8,6))
for file,label in zip(files,["GLISSER on opio (normal hist)", "MERCURY on avignon"]):
    data = np.loadtxt(file)
    t, ntp = data[:,0]/365.25, data[:,1]
    ax.plot(t, ntp, label=label)

ax.grid(linestyle="dashed")
# ax.set_xscale('log')
# ax.set_xlim(1e9,10957500000)
ax.legend()
plt.savefig("discard_6_onlyN.jpg",dpi=150)