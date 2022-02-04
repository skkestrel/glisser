import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import os
import matplotlib
import orbitplotkit as opk

font = {'weight' : 'bold',
        'size'   : 16}
matplotlib.rc('font', **font)



folder = "/home/yhuang/GLISSER/glisser/benchmark/rogue-particles-JSUNR-80000-V100/reoutput/"
output_folder = folder + "pics/hist/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)



def plotHistRes(target_pa, target_pl, k, j):
    for pa in target_pa:
        pl_txt = "hist/pl_{0:d}.txt".format(target_pl)
        pa_txt = "hist/pa_{0:d}.txt".format(pa)
        data = np.loadtxt(folder+pl_txt, usecols=[0,1,2,3,4,5,6])
        t_pl, a_pl, e_pl, I_pl, O_pl, o_pl, M_pl = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
        t_pl /= 365.25*1e6
        q_pl = a_pl * (1 - e_pl)

        data = np.loadtxt(folder+pa_txt, usecols=[0,1,2,3,4,5,6])
        t_pa, a_pa, e_pa, I_pa, O_pa, o_pa, M_pa = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
        t_pa /= 365.25*1e6
        q_pa = a_pa * (1 - e_pa)
        Q_pa = a_pa * (1 + e_pa)


        n_pl = t_pl.size
        n_pa = t_pa.size

        phi = opk.resonantAngleOuter(O_pl[0:n_pa], o_pl[0:n_pa], M_pl[0:n_pa], O_pa, o_pa, M_pa, j, k)
        loc = opk.resonanceLocationBYA(np.mean(a_pl), k, j)
        fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(5, figsize=(16,12))
        
        par_size = 2
        alpha = 1
        
        ax1.scatter(t_pa, a_pa, alpha=alpha, s=par_size)
        ax2.scatter(t_pa, q_pa, alpha=alpha, s=par_size)
        axI = ax2.twinx()
        axI.scatter(t_pa, np.rad2deg(I_pa), alpha=alpha, s=par_size,color='C1')

        ax3.scatter(t_pa, opk.wrapTo360(np.rad2deg(phi)), alpha=alpha, s=par_size)
        ax4.scatter(t_pa, opk.wrapTo360(np.rad2deg(o_pa)), alpha=alpha, s=par_size)
        axO = ax4.twinx()
        axO.scatter(t_pa, opk.wrapTo360(np.rad2deg(O_pa)), alpha=alpha, s=par_size, color='C1')

        ax5.scatter(t_pa, opk.wrapTo360(np.rad2deg(o_pa+O_pa)), alpha=alpha, s=par_size)


        # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
        # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

        aN = 30

        # ax2.plot([t_pl[0], t_pl[-1]], [aN, aN], color='red', linestyle='dashed',lw=2)
        ax2.plot([t_pl[0], t_pl[-1]], [38, 38], color='red', alpha=0.5, linestyle='dashed',lw=2)
        ax1.plot([t_pl[0], t_pl[-1]], [loc, loc], color='orange', linestyle='dashed',lw=2)
        ax4.plot([t_pl[0], t_pl[-1]], [90, 90], color='orange', linestyle='dashed',lw=2)
        ax4.plot([t_pl[0], t_pl[-1]], [270, 270], color='orange', linestyle='dashed',lw=2)


        ax1.set_title("Particle: {idx:d}".format(idx=pa))
        ax1.set_ylabel("a (au)")    
        ax2.set_ylabel("q (au)")  
        axI.set_ylabel("I (deg)")
        ax3.set_ylabel(r"$\varphi_{{{0:d}/{1:d}}}$ (deg)".format(j,k)) 
        ax3.set_ylim([0, 360])
        ax4.set_ylabel(r"$\omega$ (deg)") 
        axO.set_ylabel(r"$\Omega$ (deg)") 
        ax4.set_ylim([0, 360])
        axO.set_ylim([0, 360])

        ax5.set_ylabel(r"$\varpi$ (deg)")
        ax5.set_ylim([0, 360])
        ax5.set_xlabel("Time (Myr)")
        for ax in [ax1, ax2, ax3, ax4, ax5]:
            ax.set_xlim([t_pl[0], t_pl[-1]])


        plt.savefig(output_folder+"particle_{idx:d}_{j:d}{k:d}.jpg".format(idx=pa,j=j,k=k),dpi=200)
        print("Saved! particle: {idx:04d}".format(idx=pa))
        plt.close()
    return None

target_pl = 4
# target_pa = [44079,24744,35363,69725]
# k = 2
# j = 5

# target_pa = [29398,46339,64529,17611]
# k = 1
# j = 6

target_pa = [21956,10708,32599,34078]
k = 2
j = 9

plotHistRes(target_pa, target_pl, k, j)