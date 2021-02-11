import numpy as np
import matplotlib.pyplot as plt
import struct

folder = "/home/yhuang/GLISSER/glisser/"
filename1 = "benchmark/circular-particles-out-1000-temp/tracks/track.0.out"
# filename2 = "benchmark/circular-particles-out-ecc-read/tracks/track.0.out"

number = 20000 * 4 + 4
idx = 1
planet_list = ["Jupiter", "Saturn", "Uranus", "Neptune"]

data = np.loadtxt("temp_log_helio_1000.txt")
t, x1, y1, z1, vx1, vy1, vz1 = data[idx:number:4,0], data[idx:number:4,1], data[idx:number:4,2], data[idx:number:4,3], data[idx:number:4,4], data[idx:number:4,5], data[idx:number:4,6]
data = np.loadtxt("temp_log_jacobi_1000.txt")
t, x2, y2, z2, vx2, vy2, vz2 = data[idx:number:4,0], data[idx:number:4,1], data[idx:number:4,2], data[idx:number:4,3], data[idx:number:4,4], data[idx:number:4,5], data[idx:number:4,6]

print(t.shape)


fig, [ax1, ax2] = plt.subplots(2, figsize=(9,6))


# ax1.scatter(t, x1, s=0.2)
# ax2.scatter(t, y1, s=0.2)
# ax3.scatter(t, z1, s=0.2)
# ax4.scatter(t, vx1, s=0.2)
# ax5.scatter(t, vy1, s=0.2)
# for ax in [ax1, ax2, ax3, ax4, ax5]:
#     ax.set_xlim(0,t[-1])

for filename,label in zip([filename1],["No Hist"]):
    t, a1, e1, I1, O1, o1, F1= [], [], [],[], [], [], []
    with open(folder + filename, 'rb') as f:
        read = f.read(16)
        while len(read) == 16:
            time, npl = struct.unpack('=dQ', read)
            t.append(time)
            # print(time, npl)
            for i in range(npl):
                pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
                if(pid==idx+1):
                    a1.append(a)
                    e1.append(e)
                    I1.append(I)
                    O1.append(O)
                    o1.append(o)
                    F1.append(F)
            npa, = struct.unpack('=Q', f.read(8))
            # print(npa)
            for i in range(npa):
                pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
                # if(pid==idx):
                #     a1.append(a)
                #     e1.append(e)
                #     I1.append(I)
                #     O1.append(O)
                #     o1.append(o)
            read = f.read(16)
    print(np.array(t).shape)
    dx, dy, dz, dvx, dvy, dvz = np.array(a1) - x1, np.array(e1) - y1, np.array(I1) - z1, np.array(O1) - vx1, np.array(o1) - vy1, np.array(F1) - vz1
    dr = np.sqrt(dx**2+dy**2+dz**2)
    dv = np.sqrt(dvx**2+dvy**2+dvz**2)

    ax1.scatter(np.array(t), dr/np.sqrt(x1**2+y1**2+z1**2), alpha=0.5, label='Helio', s=0.8)
    ax2.scatter(np.array(t), dv/np.sqrt(vx1**2+vy1**2+vz1**2), alpha=0.5, label='Helio', s=0.8)

    dx, dy, dz, dvx, dvy, dvz = np.array(a1) - x2, np.array(e1) - y2, np.array(I1) - z2, np.array(O1) - vx2, np.array(o1) - vy2, np.array(F1) - vz2
    dr = np.sqrt(dx**2+dy**2+dz**2)
    dv = np.sqrt(dvx**2+dvy**2+dvz**2)

    ax1.scatter(np.array(t), dr/np.sqrt(x1**2+y1**2+z1**2), alpha=0.5, label='Jacobi', s=0.8)
    ax2.scatter(np.array(t), dv/np.sqrt(vx1**2+vy1**2+vz1**2), alpha=0.5, label='Jacobi', s=0.8)


# idarray = np.argsort(np.abs(np.diff(a1)))[::-1]
# for id in idarray[:10]:
#     print(id, t[id], np.abs(np.diff(a1))[id])
ax2.set_xlabel("Time (Days)")    

ax1.set_ylabel("dr")  
ax2.set_ylabel("dv")   

minor_ticks = np.arange(0, number * 100 +1, 100000)
major_ticks = np.arange(0, number * 100 +1, 1000000)
# minor_ticks = np.arange(0, 101, 5)


for ax in [ax1, ax2]:
    ax.set_xticks(minor_ticks,minor=True)
    ax.set_xticks(major_ticks)
    ax.set_xlim(-1000,t[-1])
    ax.set_ylim(0,0.01)
    ax.grid(which='both',alpha=0.5,linestyle='-')
    ax.legend()



plt.savefig("JACOBI_{planet}_1000.png".format(planet=planet_list[idx]),dpi=300)
# plt.show()