import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import struct
import os
import matplotlib
import orbitplotkit as opk

font = {'weight' : 'bold',
        'size'   : 16}
matplotlib.rc('font', **font)


folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Giant-Synthetic/"
target_folder = folder + "reoutput/hist/"
output_folder = folder + "pics/hist/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

target_pl = [1,2,3,4,5]


for pl in target_pl:
    pl_txt = "pl_{0:d}.txt".format(pl)
    data = np.loadtxt(target_folder+pl_txt, usecols=[0,1,2,3,4,5,6])
    t_pl, a_pl, e_pl, I_pl, O_pl, o_pl, M_pl = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
    t_pl /= 365.25*1e6
    q_pl = a_pl * (1 - e_pl)
    n_pl = t_pl.size

    fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(5, figsize=(16,12))
    
    par_size = 2
    alpha = 1
    
    ax1.scatter(t_pl, a_pl, alpha=alpha, s=par_size)
    ax2.scatter(t_pl, q_pl, alpha=alpha, s=par_size)
    axI = ax2.twinx()
    axI.scatter(t_pl, np.rad2deg(I_pl), alpha=alpha, s=par_size,color='C1')

    ax4.scatter(t_pl, opk.wrapTo360(np.rad2deg(o_pl)), alpha=alpha, s=par_size)
    axO = ax3
    axO.scatter(t_pl, opk.wrapTo360(np.rad2deg(O_pl)), alpha=alpha, s=par_size)


    node1 = (a_pl * (1 - e_pl**2))/(1+e_pl*np.cos(o_pl))
    node2 = (a_pl * (1 - e_pl**2))/(1-e_pl*np.cos(o_pl))
    ax5.scatter(t_pl, node1, alpha=alpha, s=par_size)
    ax5.scatter(t_pl, node2, alpha=alpha, s=par_size)


    # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
    # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

    aN = 30

    # ax2.plot([t_pl[0], t_pl[-1]], [aN, aN], color='red', linestyle='dashed',lw=2)
    ax2.plot([t_pl[0], t_pl[-1]], [38, 38], color='red', alpha=0.5, linestyle='dashed',lw=2)

    ax1.set_title("Particle: {idx:d}".format(idx=pl))
    ax1.set_ylabel("a (au)")    
    ax2.set_ylabel("q (au)")  
    axI.set_ylabel("I (deg)")

    ax3.set_ylim([0, 360])
    ax4.set_ylabel(r"$\omega$ (deg)") 
    axO.set_ylabel(r"$\Omega$ (deg)") 
    ax4.set_ylim([0, 360])
    axO.set_ylim([0, 360])

    # ax5.set_ylabel(r"$\varpi$ (deg)")
    ax5.set_ylabel("Nodal distance (au)")
    ax5.set_ylim([30, 100])
    ax5.set_xlabel("Time (Myr)")
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.set_xlim([t_pl[0], t_pl[-1]])


    plt.savefig(output_folder+"planet_{idx:d}.jpg".format(idx=pl),dpi=200)
    print("Saved! planet: {idx:04d}".format(idx=pl))
    plt.close()