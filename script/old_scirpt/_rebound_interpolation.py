import numpy as np
import matplotlib.pyplot as plt
import struct

folder1 = "/home/yhuang/GLISSER/glisser/examples/out-REBOUND-100/"
folder2 = "/home/yhuang/GLISSER/glisser/examples/out-REBOUND-5000/"
filename1 = "states_REBOUND.txt"
filename2 = "temp_log.out"
# filename2 = "benchmark/circular-particles-out-ecc-read/tracks/track.0.out"


fig, [ax1, ax2] = plt.subplots(2, figsize=(16,9))
data1 = np.loadtxt(folder2 + filename1)
data2 = np.loadtxt(folder1 + filename2)
for index, label, file in zip([0, 1, 2],["Uranus", "Neptune", "Rogue"], ["plr_3.txt", "plr_4.txt", "plr_5.txt"]):


    t1, x1, y1, z1 = data1[index::3,0], data1[index::3,1], data1[index::3,2], data1[index::3,3]
    t2, x2, y2, z2 = data2[index::3,0], data2[index::3,1], data2[index::3,2], data2[index::3,3]
    t = t2

    num = t.shape[0]

    data = np.loadtxt(folder2 + file, usecols = [1,5])
    a, Ome = data[:num, 0], data[:num, 1],

    print(np.array(t).shape)
    dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
    dr = np.sqrt(dx**2+dy**2+dz**2)

    ax1.scatter(t, dr, alpha=0.5, s=0.8, label = label)
    ax2.scatter(t, a,  alpha=0.5, s=0.8, label = label)
    # ax3.scatter(t, np.rad2deg(Ome),  alpha=0.5, s=0.8, label = label)

ax1.set_title("Interpolation Interval: 100 time steps")

ax2.set_xlabel("Time (Days)")    

ax1.set_ylabel("dr (au)")   
ax2.set_ylabel("a")   
# minor_ticks = np.arange(0, 1000 * 100 +1, 100000)
major_ticks = np.arange(0, int(t[-1]) +1, 100 * 1000)
# minor_ticks = np.arange(0, 101, 5)


ax1.set_xticks(major_ticks)
ax1.set_yscale('log')
ax1.set_xlim(-1e5,t[-1]/2+1e5)
ax1.set_ylim(1e-3,200)
ax1.grid(which='both',alpha=0.5,linestyle='-')


ax2.set_xticks(major_ticks)
ax2.grid(which='both',alpha=0.5,linestyle='-')
ax2.set_xlim(-1e5,t[-1]/2+1e5)
ax2.set_ylim(18,80)
ax1.legend()

plt.savefig("REBOUND_JUMP_100.png".format(index),dpi=300)
# plt.show()