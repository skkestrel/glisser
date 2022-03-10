import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import os
import matplotlib
import orbitplotkit as opk


font = {'weight': 'bold',
        'size': 16}
matplotlib.rc('font', **font)


folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNR-Resonant-4Gy-filter/"
output_folder = folder + "pics/hist/"
enc_folder = folder + "reoutput/enc/"
snapshot_folder = folder + "reoutput/snapshots/"
label = "GLISSER"
encounter_file = "encounter.out"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

if not os.path.exists(enc_folder):
    os.makedirs(enc_folder)


def resonantRatio(a_pl, a_pa):
    klist, jlist, a_rlist = opk.genOuterRes(
        a_pl, 49, 126, high1=5, high2=50, order_lim=45)
    # print(np.sort(a_rlist))

    idx = np.argmin(np.abs(a_pa-a_rlist))

    return klist[idx], jlist[idx], a_rlist[idx]


def resonantRatioRunningWindow(a_pl, a_pa, window_size=0.12, reslim=0.08, scalim=0.0375):
    code, k, j, a_res = -1, -1, -1, -1
    a_plmean, a_pamean = np.mean(a_pl), np.mean(a_pa)
    num_pl, num_pa = a_pl.shape[0], a_pa.shape[0]
    if num_pl < num_pa:
        print("Wrong planet hist size!!")
        return -1, -1, -1, [], [], []
    if num_pa < num_pl:
        print("Particle hist size is small!! Return unstable code!!")
        return -1, -1, -1, [], [], []
    if np.max(a_pa) - np.min(a_pa) > a_pamean * scalim:
        # output scattering objects:
        return -1, -1, -1, [], [], []

    num_window = int(num_pl * window_size)
    num_step = int(num_pl - num_window + 1)

    klist, jlist, a_rlist = opk.genOuterRes(
        a_plmean, 49, 126, high1=5, high2=50, order_lim=45)
    # print(np.sort(a_rlist))

    apl_window, apa_window, window_idx = [], [], []
    counter = 0
    for left in np.arange(num_step):
        right = int(left + num_window)
        window_idx.append(int(left + num_window/2))
        apl, apa = a_pl[left:right], a_pa[left:right]
        aplmean, apamean = np.mean(apl), np.mean(apa)
        if left == 0:
            k, j, a_res = resonantRatio(aplmean, apamean)
            # print(k, j, a_res)
        a_res = opk.resonanceLocationBYA(aplmean, k, j)
        if abs(apamean - a_res) < reslim:
            counter += 1

        apl_window.append(aplmean)
        apa_window.append(apamean)

    if counter/num_step > 0.5:
        # output resonant objects:
        return 1, k, j, apl_window, apa_window, window_idx
    else:
        # output non-resonant, non-scattering objects:
        k, j, a_res = resonantRatio(a_plmean, a_pamean)
        return 0, k, j, apl_window, apa_window, window_idx

        # print(num_pl, num_pa, num_window, counter/num_step ,a_plmean, a_pamean)

        # genOuterRes(a_N, 49, 126)


def plotHistRes(target_pa, target_pl, code, k, j, outputEnc):

    pl_txt = "reoutput/hist/pl_{0:d}.txt".format(target_pl)
    pa_txt = "reoutput/hist/pa_{0:d}.txt".format(target_pa)

    data = np.loadtxt(
        folder+pl_txt, usecols=[0, 1, 2, 3, 4, 5, 6])
    t_pl, a_pl, e_pl, I_pl, O_pl, o_pl, M_pl = data[:, 0], data[:,
                                                                1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6]
    t_pl /= 365.25*1e6
    q_pl = a_pl * (1 - e_pl)
    t_window = t_pl[window_idx]

    data = np.loadtxt(
        folder+pa_txt, usecols=[0, 1, 2, 3, 4, 5, 6])
    t_pa, a_pa, e_pa, I_pa, O_pa, o_pa, M_pa = data[:, 0], data[:,
                                                                1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6]
    t_pa /= 365.25*1e6
    q_pa = a_pa * (1 - e_pa)
    Q_pa = a_pa * (1 + e_pa)

    if outputEnc:
        command = "cat {2}{0} | grep ' {1} ' > {2}reoutput/enc/pa_{1}.txt".format(
            encounter_file, target_pa, folder)
        os.system(command)

        enc_txt = "reoutput/enc/pa_{0:d}.txt".format(target_pa)
        data = np.loadtxt(folder+enc_txt, usecols=[3, 5, 7])
        if len(data.shape) == 1:
            if data.shape[0] == 0:
                enc_pl, enc_t, enc_dist = [], [], []
            else:
                enc_pl, enc_t, enc_dist = [data[0]], [
                    data[1]/365.25/1e6], [data[2]]
        else:
            enc_pl, enc_t, enc_dist = data[:,
                                           0], data[:, 1]/365.25/1e6, data[:, 2]
        idxN = (enc_pl == 4)
        idxR = (enc_pl == 5)

    n_pl = t_pl.size
    n_pa = t_pa.size

    fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(5, figsize=(16, 12))

    par_size = 2
    alpha = 1

    ax1.scatter(t_pa, a_pa, alpha=0.5, s=par_size)

    if len(t_window != 0):
        ax1.scatter(t_window, apa_window, alpha=alpha, s=par_size, color='red')

    if code >= 0:
        loc = opk.resonanceLocationBYA(np.mean(a_pl), k, j)
        ax1.plot([t_pl[0], t_pl[-1]], [loc, loc],
                 color='orange', linestyle='dashed', lw=2)

    axE = ax1.twinx()
    if outputEnc and data.shape[0] != 0:
        axE.scatter(enc_t[idxR], enc_dist[idxR], marker='x', s=50, color='red')
        axE.scatter(enc_t[idxN], enc_dist[idxN],
                    marker='x', s=50, color='orange')

    ax2.scatter(t_pa, q_pa, alpha=alpha, s=par_size)
    axI = ax2.twinx()
    axI.scatter(t_pa, np.rad2deg(I_pa),
                alpha=alpha, s=par_size, color='C1')

    if code >= 0:
        phi = opk.resonantAngleOuter(
            O_pl[0:n_pa], o_pl[0:n_pa], M_pl[0:n_pa], O_pa, o_pa, M_pa, j, k)
        ax3.scatter(t_pa, opk.wrapTo360(np.rad2deg(phi)),
                    alpha=0.5, s=par_size)
        label = "{0}:{1} Resonant".format(j, k)
        textcolor = "orange"
        subsubfolder = "{0}_{1}".format(j, k)
    if code == 0:
        label = "Detached"
        textcolor = 'C0'
        subsubfolder = "{0}".format(code)
    if code < 0:
        label = "Scattering"
        textcolor = 'red'
        subsubfolder = "{0}".format(code)
    ax3.text(t_pa[-1]/2, 180, label, color=textcolor,
             ha='center', va='center', fontsize=75, weight='bold')

    ax4.scatter(t_pa, opk.wrapTo360(np.rad2deg(o_pa)),
                alpha=alpha, s=par_size)
    axO = ax4.twinx()
    axO.scatter(t_pa, opk.wrapTo360(np.rad2deg(O_pa)),
                alpha=alpha, s=par_size, color='C1')

    ax5.scatter(t_pa, opk.wrapTo360(np.rad2deg(o_pa+O_pa)),
                alpha=alpha, s=par_size)

    # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
    # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

    aN = 30

    # ax2.plot([t_pl[0], t_pl[-1]], [aN, aN], color='red', linestyle='dashed',lw=2)
    ax2.plot([t_pl[0], t_pl[-1]], [38, 38], color='red',
             alpha=0.5, linestyle='dashed', lw=2)
    ax4.plot([t_pl[0], t_pl[-1]], [90, 90],
             color='orange', linestyle='dashed', lw=2)
    ax4.plot([t_pl[0], t_pl[-1]], [270, 270],
             color='orange', linestyle='dashed', lw=2)

    ax1.set_title("Particle: {idx:d}".format(idx=pa), fontsize=30)
    ax1.set_ylabel("a (au)")
    axE.set_ylabel("Enc dist (au)")
    axE.set_ylim([0.01, 5])
    axE.set_yscale('log')
    ax2.set_ylabel("q (au)")
    axI.set_ylabel("I (deg)")
    ax3.set_ylabel(r"$\varphi_{{{0:d}/{1:d}}}$ (deg)".format(j, k))
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

    plt.tight_layout()

    subfolder = output_folder + subsubfolder
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)
    plt.savefig(subfolder+"/particle_{idx:d}.jpg".format(idx=pa), dpi=200)
    print("Saved! particle: {idx:04d}\n".format(idx=pa))

    plt.close()

    return None


target_pl = 4
output = open(folder + "classes.txt", "w")

for target_pa in np.arange(1, 100001):

    pl, pa = target_pl, target_pa
    pl_txt, pa_txt = "reoutput/hist/pl_{0:d}.txt".format(
        pl), "reoutput/hist/pa_{0:d}.txt".format(pa)
    try:
        a_pl, a_pa = np.loadtxt(
            folder+pl_txt, usecols=[1]), np.loadtxt(folder+pa_txt, usecols=[1])

        code, k, j, apl_window, apa_window, window_idx = resonantRatioRunningWindow(
            a_pl, a_pa)
        # plotHistRes(pa, pl, code, k, j, False)
        print(code, k, j, pa)
        output.write("{0:d} {1:d} {2:d} {3:d}\n".format(
            pa, code, k, j))
    except:
        print("File doesn't exist!!")
# output.close()
