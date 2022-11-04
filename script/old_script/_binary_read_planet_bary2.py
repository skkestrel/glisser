import numpy as np
import matplotlib.pyplot as plt
import struct

folder = "/home/yhuang/GLISSER/glisser/"
filename1 = "benchmark/circular-particles-out-ecc-noread/tracks/track.0.out"
filename2 = "benchmark/circular-particles-out-ecc-read/tracks/track.0.out"

files = [filename1,filename2]

# data = np.loadtxt("temp_log.txt")
# t, a, e, i, Ome, ome, M = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]



fig, [ax1, ax2, ax3, ax4, ax5, ax6] = plt.subplots(6, figsize=(9,12))
fig2, [ax7, ax8] = plt.subplots(2, figsize=(9,6))
# ax1.scatter(t, a, s=0.2)
# ax2.scatter(t, e, s=0.2)
# ax3.scatter(t, i, s=0.2)
# ax4.scatter(t, Ome, s=0.2)
# ax5.scatter(t, ome, s=0.2)
# for ax in [ax1, ax2, ax3, ax4, ax5]:
#     ax.set_xlim(0,t[-1])

idx = 1
planet_list = ["Jupiter", "Saturn", "Uranus", "Neptune"]

for filename,label in zip(files,["No Hist"]):
    t, a1, e1, I1, O1, o1, M1 = [], [], [],[], [], [], []
    with open(folder + filename, 'rb') as f:
        read = f.read(16)
        while len(read) == 16:
            time, npl = struct.unpack('=dQ', read)
            t.append(time)
            # print(time, npl)
            for i in range(npl):
                pid, a, e, I, O, o, M = struct.unpack('=I6f', f.read(28))
                if(pid==idx):
                    a1.append(a)
                    e1.append(e)
                    I1.append(I)
                    O1.append(O)
                    o1.append(o)
                    M1.append(M)
            npa, = struct.unpack('=Q', f.read(8))
            # print(npa)
            for i in range(npa):
                pid, a, e, I, O, o, M = struct.unpack('=I6f', f.read(28))
                # if(pid==idx):
                #     a1.append(a)
                #     e1.append(e)
                #     I1.append(I)
                #     O1.append(O)
                #     o1.append(o)
                #     M1..append(M)
            read = f.read(16)
    ax1.scatter(np.array(t), a1, alpha=0.5, label=label, s=0.2)
    ax2.scatter(np.array(t), e1, alpha=0.5, label=label, s=0.2)
    ax3.scatter(np.array(t), I1, alpha=0.5, label=label, s=0.2)
    ax4.scatter(np.array(t), O1, alpha=0.5, label=label, s=0.2)
    ax5.scatter(np.array(t), o1, alpha=0.5, label=label, s=0.2)
    ax6.scatter(np.array(t), M1, alpha=0.5, label=label, s=0.2)

    ax7.scatter(np.array(e1)*np.cos(np.array(o1)+np.array(O1)), np.array(e1)*np.sin(np.array(o1)+np.array(O1)), alpha=0.5, label=label, s=0.2)
    ax8.scatter(np.array(I1)*np.cos(np.array(O1)), np.array(I1)*np.sin(np.array(O1)), alpha=0.5, label=label, s=0.2)



number = len(t)
minor_ticks = np.arange(0, number * 100 +1, 1000 * 100)
major_ticks = np.arange(0, number * 100 +1, 1000 * 100 * 10)

ax5.set_xlabel("Time (Days)")    

ax1.set_ylabel("a")  
ax2.set_ylabel("e")   
ax3.set_ylabel("I")  
ax4.set_ylabel("Ome")
ax5.set_ylabel("ome")   
ax6.set_ylabel("M")   
ax1.legend()

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.set_xticks(minor_ticks,minor=True)
    ax.set_xticks(major_ticks)
    # ax.set_xlim(-1000,t[-1])
    # ax.set_ylim(0,0.01)
    ax.grid(which='both',alpha=0.5,linestyle='-')
    # ax.legend()

# for ax in [ax1, ax2, ax3, ax4, ax5]:
#     ax.set_xlim(0,1.5e8)

# ax2.legend()
planet = planet_list[idx - 1]

fig.savefig("Helio_Elements_{planet}_{id}.png".format(planet=planet, id=idx),dpi=300)
fig2.savefig("Helio_Polar_{planet}_{id}.png".format(planet=planet, id=idx),dpi=300)
# plt.show()