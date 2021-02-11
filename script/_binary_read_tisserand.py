import numpy as np
import matplotlib.pyplot as plt
import struct

folder = "/home/yhuang/GLISSER/glisser/"
filename1 = "benchmark/tisserand-particles-out-read-1/tracks/track.0.out"
filename2 = "benchmark/tisserand-particles-out-read-1000/tracks/track.0.out"

files = [filename1,filename2]
idx = 863
fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(5, figsize=(9,12))
fig2, ax6 = plt.subplots(figsize=(9,6))
for filename,label in zip(files,["1 step/block", "1000 step/block"]):
    t, a1, e1, I1, O1, o1= [], [], [], [], [], []
    with open(folder + filename, 'rb') as f:
        read = f.read(16)
        while len(read) == 16:
            time, npl = struct.unpack('=dQ', read)
            t.append(time)
            # print(time, npl)
            for i in range(npl):
                pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
                # if(pid==idx):
                #     idx = pid
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
                    idx = pid
                    a1.append(a)
                    e1.append(e)
                    I1.append(I)
                    O1.append(O)
                    o1.append(o)
            read = f.read(16)
    num = len(a1)
    print(len(t))
    t = np.array(t)[0:num]
    a = np.array(a1)
    e = np.array(e1)
    I = np.array(I1)
    size = 2
    ax1.scatter(t, a, alpha=0.5, label=label, s=size)
    ax2.scatter(t, e, alpha=0.5, label=label, s=size)
    ax3.scatter(t, np.rad2deg(I), alpha=0.5, label=label, s=size)
    ax4.scatter(t, np.rad2deg(O1), alpha=0.5, label=label, s=size)
    ax5.scatter(t, np.rad2deg(o1), alpha=0.5, label=label, s=size)

    Cj = 30/a + 2 * np.sqrt(a/30*(1-e**2))*np.cos(I)

    ax6.scatter(t, np.abs((Cj-Cj[0])/Cj[0]), alpha=0.5, label=label, s=size)

    

# idarray = np.argsort(np.abs(np.diff(a1)))[::-1]
# for id in idarray[:10]:
#     print(id, t[id], np.abs(np.diff(a1))[id])
ax5.set_xlabel("Time (days)")    
# ax1.set_xlim(5.74e8,5.75e8)
# ax2.set_xlim(5.74e8,5.75e8)
ax1.set_ylabel("a (au)")  
ax2.set_ylabel("e")    
ax3.set_ylabel("I (deg)")    
ax4.set_ylabel("Omega (deg)")    
ax5.set_ylabel("omega (deg)")    
ax1.legend()
ax6.set_yscale("log")
ax6.set_ylim(1e-12,1e-3)
# ax2.legend()
fig.savefig("tisserand3_particle_{id}_long_comp2.png".format(id=idx),dpi=200)
fig2.savefig("tisserand3_particle_{id}_Cj_comp2.png".format(id=idx),dpi=200)