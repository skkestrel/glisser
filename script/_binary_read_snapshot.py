import numpy as np
import matplotlib.pyplot as plt
import struct

folder = "/home/yhuang/GLISSER/glisser/"
filename1 = "benchmark/rogue-particles-2-out-read/tracks/track.0.out"
# filename2 = "benchmark/kuiper-particles-out-1000-read/tracks/track.0.out"

files = [filename1]
                                    
filename = filename1
label = "GLISSER"


t, a_pl, e_pl, I_pl, O_pl, o_pl= [], [], [], [], [], []
a_pa, e_pa, I_pa, O_pa, o_pa= [], [], [], [], []
with open(folder + filename, 'rb') as f:
    read = f.read(16)
    while len(read) == 16:
        time, npl = struct.unpack('=dQ', read)
        t.append(time)
        # print(time, npl)
        a1, e1, I1, O1, o1= [], [], [], [], []
        for i in range(npl):
            pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
            a1.append(a)
            e1.append(e)
            I1.append(I)
            O1.append(O)
            o1.append(o)
        a_pl.append(np.array(a1))
        e_pl.append(np.array(e1))
        I_pl.append(np.array(I1))
        O_pl.append(np.array(O1))
        o_pl.append(np.array(o1))

        npa, = struct.unpack('=Q', f.read(8))
        # print(npa)
        a1, e1, I1, O1, o1= [], [], [], [], []
        for i in range(npa):
            pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))
            a1.append(a)
            e1.append(e)
            I1.append(I)
            O1.append(O)
            o1.append(o)
        a_pa.append(np.array(a1))
        e_pa.append(np.array(e1))
        I_pa.append(np.array(I1))
        O_pa.append(np.array(O1))
        o_pa.append(np.array(o1))
        read = f.read(16)

print(len(t))

for idx in np.arange(len(t)):
    fig, [ax1, ax2] = plt.subplots(2, figsize=(9,10))
    q_pl = a_pl[idx] * (1 - e_pl[idx])
    ax1.scatter(a_pl[idx], q_pl, alpha=0.8, label=label, s=100, color='red', marker='x',lw=2)
    ax2.scatter(a_pl[idx], np.rad2deg(I_pl[idx]), alpha=0.8, label=label, s=100, color='red', marker='x',lw=2)

    q_pa = a_pa[idx] * (1 - e_pa[idx])
    ax1.scatter(a_pa[idx], q_pa, alpha=0.5, label=label, s=5)
    ax2.scatter(a_pa[idx], np.rad2deg(I_pa[idx]), alpha=0.5, label=label, s=5)

    outer = 600
    x = np.linspace(5,outer,1000)
    ax1.plot(x,x, color='grey', linestyle='dashed',lw=1)


    ax2.set_xlabel("a (au)")    
    ax1.set_ylabel("q (au)")  
    ax2.set_ylabel("I (deg)")    

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax1.set_xlim(25,700)
    ax1.set_ylim(25,200)
    ax2.set_xlim(25,700)
    ax2.set_ylim(0,60)

    plt.title("Time: {time:6.3f} Myr".format(time=t[idx]/365.25/1e6))

    plt.savefig("Rogue2/snapshot_{idx:04d}.jpg".format(idx=idx),dpi=120)
    print("Saved! frame: {idx:04d}".format(idx=idx))