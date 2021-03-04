import numpy as np
import matplotlib.pyplot as plt
import struct

folder = "/home/yhuang/GLISSER/glisser/"
filename1 = "benchmark/rogue-particles-30m-onlyN/tracks/track.0.out"
# filename1 = "sockeye_results/track_100.out"
# filename2 = "benchmark/kuiper-particles-out-1000-read/tracks/track.0.out"

files = [filename1]
                                    
filename = filename1
label = "GLISSER"


t, a_pl, e_pl, I_pl, O_pl, o_pl= [], [], [], [], [], []
a_pa, e_pa, I_pa, O_pa, o_pa= [], [], [], [], []
file = open("discard_opio_onlyN_5.out", "w+")
with open(folder + filename, 'rb') as f:
    read = f.read(16)
    while len(read) == 16:
        time, npl = struct.unpack('=dQ', read)
        t.append(time)
        # print(time, npl)
        a1, e1, I1, O1, o1= [], [], [], [], []
        for i in range(npl):
            pid, a, e, I, O, o, F = struct.unpack('=I6d', f.read(52))
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
            pid, a, e, I, O, o, F = struct.unpack('=I6d', f.read(52))
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
        file.write(str(int(time)))
        file.write(' ')
        file.write(str(int(1000-npa)))
        file.write('\n')


print(len(t))


# for idx in np.arange(len(t),step=50):
#     n_bins = 200
#     fig, [ax1, ax2, ax3] = plt.subplots(3, figsize=(14,10))
#     # q_pl = a_pl[idx] * (1 - e_pl[idx])
#     # ax1.scatter(a_pl[idx], q_pl, alpha=0.8, label=label, s=100, color='red', marker='x',lw=2)
#     # ax2.scatter(a_pl[idx], I_pl[idx], alpha=0.8, label=label, s=100, color='red', marker='x',lw=2)

#     q_pa = a_pa[idx] * (1 - e_pa[idx])
#     ax1.hist(a_pa[idx], n_bins, density=True, histtype='step',cumulative=True, label='a')
#     ax2.hist(q_pa, n_bins, density=True, histtype='step',cumulative=True, label='q')
#     ax3.hist(np.rad2deg(I_pa[idx]), n_bins, density=True, histtype='step',cumulative=True, label='I')


#     ax1.set_xlabel("a (au)")  
#     ax1.set_xlim(25,2000)
#     ax1.set_xscale("log")
#     ax2.set_xlabel("q (au)")
#     ax2.set_xlim(10,120)
#     ax3.set_xlabel("I (deg)")
#     ax3.set_xlim(0,40)

#     ax1.grid(which="major", linestyle="dashed")
#     ax1.grid(which="minor", linestyle="dotted")
#     ax2.grid(linestyle="dashed")
#     ax3.grid(linestyle="dashed")

#     ax1.set_title("Time: {time:6.3f} Myr".format(time=t[idx]/365.25/1e6))

#     plt.savefig("Rogue_histograms4/histo_{idx:04d}.jpg".format(idx=idx),dpi=120)
#     print("Saved! frame: {idx:04d}".format(idx=idx))