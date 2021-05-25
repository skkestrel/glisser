import numpy as np
import matplotlib.pyplot as plt
import struct

filename1 = "benchmark/rogue-particles-JSUNR-1000/discard.out"
filename2 = "sockeye_results/discard_1000.out"
filename3 = "avignon_results/info.out"
filename4 = "benchmark/rogue-particles-JSUNR-80000/discard.out"
filename5 = "sockeye_results/discard_8cores.out"
filename6 = "sockeye_results/discard_12cores.out"


# filename2 = "benchmark/kuiper-particles-out-1000-read/tracks/track.0.out"

files = [filename1, filename2,filename3,filename4,filename5,filename6]
# files = [filename1]

fig, ax = plt.subplots(figsize=(10,6))

for file, totalnum, label in zip(files,[1e3,1e3,1e3,8e4,8e4,8e4],["GLISSER/1k/opio","GLISSER/1k/sockeye", "MERCURY/1k/avignon", "GLISSER/80k/opio", "GLISSER/80k/sockeye/8cores", "GLISSER/80k/sockeye/12cores"]):
    if label[0] == "G":
        data = np.loadtxt(file, usecols=[1,-1])
    else:
        data = np.loadtxt(file, usecols=[0,-2])
    code, t = data[:,0], data[:,1]/365.25/1e6
    if totalnum < 1e4:
        ax.plot(t, (np.arange(t.shape[0])+1)/totalnum, label=label, alpha=0.75, lw=2)
    else:
        ax.plot(t, (np.arange(t.shape[0])+1)/totalnum, label=label, alpha=0.75, lw=2)

# ax1.hist(a_pa[idx], n_bins, density=True, histtype='step',cumulative=True, label='a')
# ax2.hist(q_pa, n_bins, density=True, histtype='step',cumulative=True, label='q')


ax.grid(linestyle="dashed")
# ax.set_xscale('log')
ax.set_xlabel("Time (Myr)")  
ax.set_ylabel("Ratio of Discarded Particles")  
ax.set_xlim(0.9,105)
ax.legend()
plt.savefig("discard_JSUNR_new.png",dpi=150)