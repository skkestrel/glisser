import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import struct
import os
import matplotlib
import orbitplotkit as opk

# font = {'weight' : 'bold',
#         'size'   : 14}
# matplotlib.rc('font', **font)


folder = "/home/yhuang/GLISSER/glisser/examples/out-JSUNR-Resonant-Test-4000-bary/reoutput/"
output_folder = folder + "pics/hist/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

target_pl = [1]
# target_pa = np.arange(3000, 4001, 10)


for pl in target_pl:

    pl_txt = "hist/pl_{0:d}.txt".format(pl)
    data = np.loadtxt(folder+pl_txt, usecols=[0,1,2,3,4,5,6])
    t_pl, a_pl, e_pl, I_pl, O_pl, o_pl, M_pl = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
    t_pl /= 365.25*1e6
    q_pl = a_pl * (1-e_pl)
    
    plr_txt = "hist/plr_{0:d}.txt".format(pl)
    data = np.loadtxt(folder+plr_txt, usecols=[0,1,2,3,4,5,6])
    t_pa, a_pa, e_pa, I_pa, O_pa, o_pa, M_pa = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
    t_pa /= 365.25*1e6
    q_pa = a_pa * (1-e_pa)

    n_pl = t_pl.size
    n_pa = t_pa.size

    
    fig, [ax1, ax2, ax3, ax4] = plt.subplots(4, figsize=(16,10))
    
    par_size = 2
    alpha = 1
    
    ax1.scatter(t_pl, a_pl, alpha=alpha, s=par_size, label="GLISSER")
    ax2.scatter(t_pl, q_pl, alpha=alpha, s=par_size)
    ax3.scatter(t_pl, np.rad2deg(I_pl), alpha=alpha, s=par_size)
    ax4.scatter(t_pl, opk.wrapTo360(np.rad2deg(O_pl)), alpha=alpha, s=par_size)


    ax1.scatter(t_pa, a_pa, alpha=alpha, s=par_size, label="REBOUND")
    ax2.scatter(t_pa, q_pa, alpha=alpha, s=par_size)
    ax3.scatter(t_pa, np.rad2deg(I_pa), alpha=alpha, s=par_size)
    ax4.scatter(t_pa, opk.wrapTo360(np.rad2deg(O_pa)), alpha=alpha, s=par_size)

    ax1.set_ylabel("a (au)")    
    ax2.set_ylabel("q (au)")  
    ax3.set_ylabel("I (deg)")
    ax4.set_ylabel("Omega (deg)") 
    ax4.set_ylim([0, 360])   
    ax4.set_xlabel("Time (Myr)")
    
    for ax in [ax1, ax2, ax3, ax4]:
        ax.grid(True)
        ax.set_xlim([t_pl[0], t_pl[-1]])

    ax1.legend()  
    plt.savefig(output_folder+"pl_{idx:d}.jpg".format(idx=pl),dpi=200)
    print("Saved! pl: {idx:04d}".format(idx=pl))
    plt.close()