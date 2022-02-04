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
folder_swift = "/home/yhuang/swift4/_test/"
output_folder = folder + "pics/hist/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

target_pl = 4
# target_pa = np.arange(3000, 4001, 10)
target_pa = [200,1180,2340,3220,3960]


for pa in target_pa:
    if pa <= 1000:
        j = 3
        k = 2
    elif pa <= 2000:
        j = 2
        k = 1
    elif pa <= 3000:
        j = 5
        k = 2
    else:
        j = 3
        k = 1

    print(j, k)
    pl_txt = "hist/pl_{0:d}.txt".format(target_pl)
    pa_txt = "hist/pa_{0:d}.txt".format(pa)
    data = np.loadtxt(folder+pl_txt, usecols=[0,1,2,3,4,5,6])
    t_pl, a_pl, e_pl, I_pl, O_pl, o_pl, M_pl = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
    t_pl /= 365.25*1e6
    

    data = np.loadtxt(folder+pa_txt, usecols=[0,1,2,3,4,5,6])
    t_pa, a_pa, e_pa, I_pa, O_pa, o_pa, M_pa = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
    t_pa /= 365.25*1e6
    q_pa = a_pa * (1-e_pa)

    n_pl = t_pl.size
    n_pa = t_pa.size

    phi = opk.resonantAngleOuter(O_pl[0:n_pa], o_pl[0:n_pa], M_pl[0:n_pa], O_pa, o_pa, M_pa, j, k)
    
    fig, [ax1, ax2, ax3, ax4] = plt.subplots(4, figsize=(16,10))
    
    par_size = 2
    alpha = 1
    
    ax1.scatter(t_pa, a_pa, alpha=alpha, s=par_size, label="GLISSER")
    ax2.scatter(t_pa, q_pa, alpha=alpha, s=par_size)
    ax3.scatter(t_pa, np.rad2deg(I_pa), alpha=alpha, s=par_size)
    ax4.scatter(t_pa, opk.wrapTo360(np.rad2deg(phi)), alpha=alpha, s=par_size)


    pl_txt = "pl_{0:d}.out".format(target_pl)
    pa_txt = "pa_{0:d}.out".format(pa)



    data = np.loadtxt(folder_swift+pl_txt, usecols=[0,1,2,3,4,5,6])
    t_pl, a_pl, e_pl, I_pl, O_pl, o_pl, M_pl = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
    t_pl /= 365.25*1e6

    data = np.loadtxt(folder_swift+pa_txt, usecols=[0,1,2,3,4,5,6])
    t_pa, a_pa, e_pa, I_pa, O_pa, o_pa, M_pa = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
    t_pa /= 365.25*1e6
    q_pa = a_pa * (1-e_pa)

    n_pl = t_pl.size
    n_pa = t_pa.size

    phi = opk.resonantAngleOuter(O_pl[0:n_pa], o_pl[0:n_pa], M_pl[0:n_pa], O_pa, o_pa, M_pa, j, k)

    ax1.scatter(t_pa, a_pa, alpha=alpha, s=par_size, label="SWIFT")
    ax2.scatter(t_pa, q_pa, alpha=alpha, s=par_size)
    ax3.scatter(t_pa, I_pa, alpha=alpha, s=par_size)
    ax4.scatter(t_pa, opk.wrapTo360(phi), alpha=alpha, s=par_size)

    aN = np.mean(a_pl)
    loc = opk.resonanceLocationByA(aN, k, j)
    ax2.plot([t_pl[0], t_pl[-1]], [aN, aN], color='C4', linestyle='dashed',lw=2)
    ax2.plot([t_pl[0], t_pl[-1]], [38, 38], color='C4', alpha=0.5, linestyle='dashed',lw=2)
    ax1.plot([t_pl[0], t_pl[-1]], [loc, loc], color='C4', linestyle='dashed',lw=2)

    ax1.set_title("Particle: {idx:d}".format(idx=pa))
    ax1.set_ylabel("a (au)")    
    ax2.set_ylabel("q (au)")  
    ax3.set_ylabel("I (deg)")
    ax4.set_ylabel("phi_{0:d}/{1:d} (deg)".format(j,k)) 
    ax4.set_ylim([0, 360])   
    ax4.set_xlabel("Time (Myr)")
    
    for ax in [ax1, ax2, ax3, ax4]:
        ax.grid(True)
        ax.set_xlim([t_pl[0], t_pl[-1]])

    ax1.legend()  
    plt.savefig(output_folder+"particle_{idx:d}.jpg".format(idx=pa),dpi=200)
    print("Saved! particle: {idx:04d}".format(idx=pa))
    plt.close()