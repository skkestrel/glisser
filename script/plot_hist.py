import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import struct
import os
import matplotlib
import orbitplotkit as opk

font = {'weight': 'bold',
        'size': 14}
matplotlib.rc('font', **font)


folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
output_folder = folder.format(5000) + "pics/hist/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

target_pl = 4
target_pa = np.arange(10, 97461, 1000)
# target_pa = [3000]


for pa in target_pa:
    fig, [ax1, ax2, ax3, ax4] = plt.subplots(4, figsize=(16, 10))

    pl_txt = "hist/pl_{0:d}.txt".format(target_pl)
    pa_txt = "hist/pa_{0:d}.txt".format(pa)
    data = np.loadtxt(folder+pl_txt,
                      usecols=[0, 1, 2, 3, 4, 5, 6])
    t_pl, a_pl, e_pl, I_pl, O_pl, o_pl, M_pl = data[:, 0], data[:,
                                                                1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6]
    t_pl /= 365.25*1e6
    q_pl = a_pl * (1 - e_pl)

    data = np.loadtxt(folder+pa_txt,
                      usecols=[0, 1, 2, 3, 4, 5, 6])
    print(len(data.shape))
    if len(data.shape) == 1:
        continue

    num = data.shape[0]
    t_pa, a_pa, e_pa, I_pa, O_pa, o_pa, M_pa = data[:, 0], data[:,
                                                                1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6]

    t_pa /= 365.25*1e6
    q_pa = a_pa * (1 - e_pa)
    Q_pa = a_pa * (1 + e_pa)

    n_pl = t_pl.size
    n_pa = t_pa.size

    par_size = 2
    alpha = 1

    ax1.plot(t_pa, a_pa, alpha=alpha, lw=par_size)
    ax1.plot(t_pa, q_pa, alpha=alpha, lw=par_size)
    ax1.plot(t_pa, Q_pa, alpha=alpha, lw=par_size)

    ax2.plot(t_pa, np.rad2deg(I_pa), alpha=alpha, lw=par_size)
    ax3.plot(t_pa, opk.wrapTo360(np.rad2deg(O_pa + o_pa - opk.g_freqs[7]*opk.ARCSEC_TO_RAD*t_pa*1e6)), alpha=alpha,
                lw=par_size)
    ax4.plot(t_pa, opk.wrapTo360(np.rad2deg(O_pa - opk.f_freqs[7]*opk.ARCSEC_TO_RAD*t_pa*1e6)), alpha=alpha,
                lw=par_size)

    # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
    # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

    ax1.set_title("Particle: {idx:d}".format(idx=pa))
    ax1.set_ylabel("a/q/Q (au)")
    ax2.set_ylabel("I (deg)")
    ax3.set_ylabel(r"$\varpi - g_8 t$")
    ax4.set_ylabel(r"$\Omega - f_8 t$")
    ax4.set_xlabel("Time (Myr)")

    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xlim(1, t_pl[-1])
        ax.set_xscale('log')
    for ax in [ax3, ax4]:
        ax.set_ylim(0, 360)

    # ax4.legend()
    if num == t_pl.shape[0]:
        plt.savefig(output_folder+"particle_{idx:d}.jpg".format(idx=pa), dpi=200)
    else:
        plt.savefig(output_folder+"particle_{idx:d}_U.jpg".format(idx=pa), dpi=200)
    print("Saved! particle: {idx:06d}".format(idx=pa))
    plt.close()
