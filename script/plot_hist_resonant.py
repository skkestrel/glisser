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


folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNR-Scattering/"
output_folder = folder + "pics/hist/"
snapshot_folder = folder + "reoutput/snapshots/"
label = "GLISSER"
encounter_file = "encounter.out"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


def plotHistRes(target_pa, target_pl, k, j):
    command = "cat {2}{0} | grep ' {1} ' > {2}reoutput/enc/pa_{1}.txt".format(
        encounter_file, target_pa, folder)
    os.system(command)

    pl_txt = "reoutput/hist/pl_{0:d}.txt".format(target_pl)
    pa_txt = "reoutput/hist/pa_{0:d}.txt".format(target_pa)
    enc_txt = "reoutput/enc/pa_{0:d}.txt".format(target_pa)
    data = np.loadtxt(
        folder+pl_txt, usecols=[0, 1, 2, 3, 4, 5, 6])
    t_pl, a_pl, e_pl, I_pl, O_pl, o_pl, M_pl = data[:, 0], data[:,
                                                                1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6]
    t_pl /= 365.25*1e6
    q_pl = a_pl * (1 - e_pl)

    data = np.loadtxt(
        folder+pa_txt, usecols=[0, 1, 2, 3, 4, 5, 6])
    t_pa, a_pa, e_pa, I_pa, O_pa, o_pa, M_pa = data[:, 0], data[:,
                                                                1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6]
    t_pa /= 365.25*1e6
    q_pa = a_pa * (1 - e_pa)
    Q_pa = a_pa * (1 + e_pa)

    data = np.loadtxt(folder+enc_txt, usecols=[3, 5, 7])
    if len(data.shape) == 1:
        if data.shape[0] == 0:
            enc_pl, enc_t, enc_dist = [], [], []
        else:
            enc_pl, enc_t, enc_dist = [data[0]], [
                data[1]/365.25/1e6], [data[2]]
    else:
        enc_pl, enc_t, enc_dist = data[:, 0], data[:, 1]/365.25/1e6, data[:, 2]
    idxN = (enc_pl == 4)
    idxR = (enc_pl == 5)

    n_pl = t_pl.size
    n_pa = t_pa.size

    phi = opk.resonantAngleOuter(
        O_pl[0:n_pa], o_pl[0:n_pa], M_pl[0:n_pa], O_pa, o_pa, M_pa, j, k)
    loc = opk.resonanceLocationBYA(np.mean(a_pl), k, j)

    fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(5, figsize=(16, 12))

    par_size = 2
    alpha = 1

    ax1.scatter(t_pa, a_pa, alpha=alpha, s=par_size)
    axE = ax1.twinx()
    if data.shape[0] != 0:
        axE.scatter(enc_t[idxR], enc_dist[idxR], marker='x', s=50, color='red')
        axE.scatter(enc_t[idxN], enc_dist[idxN],
                    marker='x', s=50, color='orange')

    ax2.scatter(t_pa, q_pa, alpha=alpha, s=par_size)
    axI = ax2.twinx()
    axI.scatter(t_pa, np.rad2deg(I_pa),
                alpha=alpha, s=par_size, color='C1')

    ax3.scatter(t_pa, opk.wrapTo360(np.rad2deg(phi)),
                alpha=alpha, s=par_size)
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
    ax1.plot([t_pl[0], t_pl[-1]], [loc, loc],
             color='orange', linestyle='dashed', lw=2)
    ax4.plot([t_pl[0], t_pl[-1]], [90, 90],
             color='orange', linestyle='dashed', lw=2)
    ax4.plot([t_pl[0], t_pl[-1]], [270, 270],
             color='orange', linestyle='dashed', lw=2)

    ax1.set_title("Particle: {idx:d}".format(idx=pa))
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

    subfolder = output_folder + "/{j:d}{k:d}/".format(j=j, k=k)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)
    plt.savefig(
        subfolder+"particle_{idx:d}.jpg".format(idx=pa), dpi=200)
    print("Saved! particle: {idx:04d}\n".format(idx=pa))

    plt.close()

    return None


final_snapshot = "planets_36520000000.txt"
df_pl = pd.read_csv(snapshot_folder+final_snapshot, names=[
    'id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ').set_index('id')

final_snapshot = "particles_36520000000.txt"
df_pa = pd.read_csv(snapshot_folder+final_snapshot, names=[
    'id', 'a', 'e', 'inc', 'Omega', 'omega', 'f'], delimiter=' ').set_index('id')
df_pa['q'] = df_pa['a']*(1-df_pa['e'])
df_q = df_pa[df_pa['q'] > 45]
idx = np.array(df_q.index)
a_N = df_pl.iloc[3]['a']


def resonantRatio(target_pa):
    a_94 = opk.resonanceLocationBYA(a_N, 4, 9)
    a_52 = opk.resonanceLocationBYA(a_N, 2, 5)
    a_72 = opk.resonanceLocationBYA(a_N, 2, 7)
    a_92 = opk.resonanceLocationBYA(a_N, 2, 9)
    a_112 = opk.resonanceLocationBYA(a_N, 2, 11)
    a_132 = opk.resonanceLocationBYA(a_N, 2, 13)
    a_73 = opk.resonanceLocationBYA(a_N, 3, 7)
    a_83 = opk.resonanceLocationBYA(a_N, 3, 8)
    a_103 = opk.resonanceLocationBYA(a_N, 3, 10)
    a_113 = opk.resonanceLocationBYA(a_N, 3, 11)
    a_133 = opk.resonanceLocationBYA(a_N, 3, 13)
    a_143 = opk.resonanceLocationBYA(a_N, 3, 14)
    a_163 = opk.resonanceLocationBYA(a_N, 3, 16)
    a_173 = opk.resonanceLocationBYA(a_N, 3, 17)
    a_193 = opk.resonanceLocationBYA(a_N, 3, 19)
    a_203 = opk.resonanceLocationBYA(a_N, 3, 20)
    a_31 = opk.resonanceLocationBYA(a_N, 1, 3)
    a_41 = opk.resonanceLocationBYA(a_N, 1, 4)
    a_51 = opk.resonanceLocationBYA(a_N, 1, 5)
    a_61 = opk.resonanceLocationBYA(a_N, 1, 6)
    a_71 = opk.resonanceLocationBYA(a_N, 1, 7)
    a_81 = opk.resonanceLocationBYA(a_N, 1, 8)
    a_152 = opk.resonanceLocationBYA(a_N, 2, 15)

    a_rlist = np.array([a_94, a_52, a_72, a_92, a_112, a_132, a_73, a_83, a_103, a_113,
                        a_133, a_143, a_163, a_173, a_193, a_203,
                        a_31, a_41, a_51, a_61, a_71, a_81, a_152])
    klist = [4, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
             3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 2]
    jlist = [9, 5, 7, 9, 11, 13, 7, 8, 10, 11, 13,
             14, 16, 17, 19, 20, 3, 4, 5, 6, 7, 8, 15]
    target_pa

    pa_txt = "reoutput/hist/pa_{0:d}.txt".format(target_pa)

    a_pl = np.loadtxt(folder+pa_txt, usecols=[1])
    a_mean = np.mean(a_pl)
    # print(np.abs(a_mean-a_rlist))

    idx = np.argmin(np.abs(a_mean-a_rlist))

    return klist[idx], jlist[idx], a_mean, a_rlist[idx]


target_pl = 4
# target_pa = [44079,24744,35363,69725]
# k = 2
# j = 5

# target_pa = [29398,46339,64529,17611]
# k = 1
# j = 6

print(idx)
target_pa = idx

for pa in target_pa:
    if pa > 62935:
        k, j, a_mean, a_r = resonantRatio(pa)
        print(pa)
        print(k, j)
        plotHistRes(pa, target_pl, k, j)
