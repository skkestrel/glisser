import numpy as np
import matplotlib.pyplot as plt
import struct

folder = "/home/yhuang/GLISSER/glisser/"
filename1 = "benchmark/rogue-particles-1000-out-read/tracks/track.0.out"
# filename2 = "benchmark/kuiper-particles-out-1000-read/tracks/track.0.out"

files = [filename1]
idx = 482                                    
fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(5, figsize=(9,12))
for filename,label in zip(files,["GLISSER"]):
    t, a1, e1, I1, O1, o1, f1 = [], [], [], [], [], [], []
    with open(folder + filename, 'rb') as f:
        read = f.read(16)
        while len(read) == 16:
            time, npl = struct.unpack('=dQ', read)
            t.append(time)
            for i in range(npl):
                pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
                # if(pid==idx):
                #     idx = pid
                #     a1.append(a)
                #     e1.append(e)
                #     I1.append(I)
                #     O1.append(O)
                #     o1.append(o)
                #     f1.append(F)
            npa, = struct.unpack('=Q', f.read(8))
            # print(npa)
            for i in range(npa):
                pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
                if(pid==idx):
                    idx = pid
                    a1.append(a)
                    e1.append(e)
                    I1.append(I)
                    O1.append(O)
                    o1.append(o)
            read = f.read(16)
    t = np.array(t)
    O1 = np.array(O1)
    o1 = np.array(o1)
    ax1.scatter(t, a1, alpha=0.5, label=label, s=0.5)
    ax2.scatter(t, e1, alpha=0.5, label=label, s=0.5)
    ax3.scatter(t, np.rad2deg(I1), alpha=0.5, label=label, s=0.5)
    ax4.scatter(t, np.rad2deg(o1), alpha=0.5, label=label, s=0.5)
    ax5.scatter(t, np.rad2deg(O1), alpha=0.5, label=label, s=0.5)

    

# idarray = np.argsort(np.abs(np.diff(a1)))[::-1]
# for id in idarray[:10]:
#     print(id, t[id], np.abs(np.diff(a1))[id])
ax5.set_xlabel("Time (days)")    
# ax1.set_xlim(5.74e8,5.75e8)
# ax2.set_xlim(5.74e8,5.75e8)
ax1.set_ylabel("a (au)")  
ax2.set_ylabel("e")    
ax3.set_ylabel("I (deg)")    
ax4.set_ylabel("omega (deg)")    
ax5.set_ylabel("True (deg)")    
ax1.legend()

# for ax in [ax1, ax2, ax3, ax4, ax5]:
    # ax.set_xlim(-100,t[-1]/48)

# ax2.legend()
plt.savefig("test_{id}_1000.png".format(id=idx),dpi=200)
# plt.show()