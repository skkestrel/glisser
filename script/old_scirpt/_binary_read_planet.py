import numpy as np
import matplotlib.pyplot as plt
import struct

folder = "/home/yhuang/GLISSER/glisser/"
filename1 = "benchmark/circular-particles-out-ecc-noread/tracks/track.0.out"
filename2 = "benchmark/circular-particles-out-ecc-read/tracks/track.0.out"

# data = np.loadtxt("temp_log.txt",max_rows=200000)
# t, a, e, i, Ome, ome, M = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]

files = [filename1, filename2]

fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(5, figsize=(9,12))

# ax1.scatter(t, a, s=0.2)
# ax2.scatter(t, e, s=0.2)
# ax3.scatter(t, i, s=0.2)
# ax4.scatter(t, Ome, s=0.2)
# ax5.scatter(t, ome, s=0.2)
# for ax in [ax1, ax2, ax3, ax4, ax5]:
#     ax.set_xlim(0,t[-1])

idx = 1

for filename,label in zip(files,["No Hist", "Use Hist"]):
    t, a1, e1, I1, O1, o1= [], [], [],[], [], []
    with open(folder + filename, 'rb') as f:
        read = f.read(16)
        while len(read) == 16:
            time, npl = struct.unpack('=dQ', read)
            t.append(time)
            # print(time, npl)
            for i in range(npl):
                pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
                # if(pid==idx):
                #     a1.append(a)
                #     e1.append(e)
                #     I1.append(I)
                #     O1.append(O)
                #     o1.append(o)
            npa, = struct.unpack('=Q', f.read(8))
            # print(npa)
            for i in range(npa):
                pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
                if(pid==idx):
                    a1.append(a)
                    e1.append(e)
                    I1.append(I)
                    O1.append(O)
                    o1.append(o)
            read = f.read(16)
    ax1.scatter(np.array(t), a1, alpha=0.5, label=label, s=0.2)
    ax2.scatter(np.array(t), e1, alpha=0.5, label=label, s=0.2)
    ax3.scatter(np.array(t), I1, alpha=0.5, label=label, s=0.2)
    ax4.scatter(np.array(t), O1, alpha=0.5, label=label, s=0.2)
    ax5.scatter(np.array(t), o1, alpha=0.5, label=label, s=0.2)


# idarray = np.argsort(np.abs(np.diff(a1)))[::-1]
# for id in idarray[:10]:
#     print(id, t[id], np.abs(np.diff(a1))[id])
ax5.set_xlabel("Time (Days)")    

ax1.set_ylabel("a")  
ax2.set_ylabel("e")   
ax3.set_ylabel("I")  
ax4.set_ylabel("O")
ax5.set_ylabel("o")   
ax1.legend()

# for ax in [ax1, ax2, ax3, ax4, ax5]:
#     ax.set_xlim(0,1.5e8)

# ax2.legend()
plt.savefig("Ecc_particle_done_{id}.png".format(id=idx),dpi=300)
# plt.show()