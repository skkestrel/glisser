import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import struct
import os
import matplotlib

font = {'weight' : 'bold',
        'size'   : 14}
matplotlib.rc('font', **font)


folder = "/home/yhuang/GLISSER/glisser/examples/out-JSUN-Polar-100m-80000/reoutput/"
output_folder = folder + "pics/snapshots/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

time = [0]
# time = [6275000000]
time = np.arange(0,6275000001,25000000)
idx = 1


for t in time:
    pl_txt = "snapshots/planets_{0:d}.txt".format(t)
    pa_txt = "snapshots/particles_{0:d}.txt".format(t)
    data = np.loadtxt(folder+pl_txt, usecols=[1,2,3])
    a_pl, e_pl, I_pl = data[:,0], data[:,1], data[:,2]
    data = np.loadtxt(folder+pa_txt, usecols=[1,2,3])
    a_pa, e_pa, I_pa = data[:,0], data[:,1], data[:,2]

    fig, [ax1, ax2] = plt.subplots(2, figsize=(9,10))
    q_pl = a_pl * (1 - e_pl)
    ax1.scatter(a_pl, q_pl, alpha=0.8, s=100, color='red', marker='x',lw=2, zorder=2)
    ax2.scatter(a_pl, np.rad2deg(I_pl), alpha=0.8, s=100, color='red', marker='x',lw=2, zorder=2)

    q_pa = a_pa * (1 - e_pa)
    par_size = 0.1
    ax1.scatter(a_pa, q_pa, alpha=0.4, s=par_size)
    ax2.scatter(a_pa, np.rad2deg(I_pa), alpha=0.4, s=par_size)
    ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
    # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

    aN = a_pl[-2]
    outer = 100
    inner = 3
    x = np.linspace(inner,outer,1000)
    y = np.linspace(inner,outer*1000,1000)
    ax1.plot(x, x, color='grey', linestyle='dashed',lw=1)
    ax1.fill_between(x, x, y, facecolor='gray', alpha=0.5)

    ax1.plot([aN, outer], [aN, aN], color='red', linestyle='dashed',lw=2)
    ax1.plot([38, outer], [38, 38], color='red', alpha=0.5, linestyle='dashed',lw=2)
    # ax1.plot([80, outer], [80, 80], color='red', alpha=0.3, linestyle='dashed',lw=2)
    # for pos,label in zip([aN,38], ["q = 30", "q = 38"]):
    #     ax1.text(15, pos-2, label, fontsize = 12, color='red')

    rect2 = patches.Rectangle((4, 80), 46, 20, edgecolor='none', facecolor='green', alpha=0.2)
    ax2.add_patch(rect2)
    # ress = [11/5, 7/3, 5/2, 8/3, 11/4, 3/1, 13/4, 7/2, 15/4, 4/1, 9/2, 5/1]
    # labels = ["11/5", "7/3", "5/2", "8/3", "11/4", "3/1", "13/4","7/2", "15/4", "4/1","9/2", "5/1"]
    # for res, label in zip(ress, labels):
    #     loc = aN*(res)**(2/3)
    #     ax1.plot([loc, loc], [0, loc], color='grey', alpha=0.4 ,linestyle='dashed',lw=1, zorder=-1)
    #     ax2.plot([loc, loc], [0, 120], color='grey', alpha=0.4 ,linestyle='dashed',lw=1, zorder=-1)
    #     ax1.text(loc-0.6, 66.5, label, fontsize = 12, rotation=60)

    ax2.set_xlabel("a (au)")    
    ax1.set_ylabel("q (au)")  
    ax2.set_ylabel("I (deg)")    

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax1.set_xlim(inner,outer)
    ax1.set_ylim(inner,350)
    ax2.set_xlim(inner,outer)
    ax2.set_ylim(0,180)

    plt.title("Time: {time:6.3f} Myr".format(time=t/365/1e6))

    plt.savefig(output_folder+"zoom_{idx:04d}.jpg".format(idx=idx),dpi=200)
    print("Saved! frame: {idx:04d}".format(idx=idx))
    plt.close()
    idx += 1 