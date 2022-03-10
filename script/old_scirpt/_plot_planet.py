import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("temp_log.txt")
t, a, e, i, Ome, ome, M = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]

M = np.fmod(M, 2*np.pi)

fig, [ax1, ax2, ax3, ax4, ax5, ax6] = plt.subplots(6, figsize=(9,12))

ax1.scatter(t, a, s=0.2)
ax2.scatter(t, e, s=0.2)
ax3.scatter(t, i, s=0.2)
ax4.scatter(t, Ome, s=0.2)
ax5.scatter(t, ome, s=0.2)
ax6.scatter(t, M, s=0.2)

# minor_ticks = np.arange(0,t[-1]+1,10000)
# print(minor_ticks)
# ax5.set_xticks(minor_ticks,minor=True)
# ax5.grid(which='minor', alpha=0.2)
# ax5.grid(which='major', alpha=0.5)
# ax6.set_xticks(minor_ticks,minor=True)
# ax6.grid(which='minor', alpha=0.2)
# ax5.grid(which='major', alpha=0.5)

# for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
#     ax.set_xlim(3.5560e8,3.563e8)


ax3.set_xlabel("Time (Days)")    
# ax.set_ylabel("energy")  
# ax.legend()
plt.savefig("Neptune_bary.png",dpi=300)