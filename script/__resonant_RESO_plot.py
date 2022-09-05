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


folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Giant-Synthetic/"
output_folder = folder + "pics/hist/"
enc_folder = folder + "reoutput/enc/"
histRESO_folder = folder + "reoutput/hist_RESO/"
hist_folder = folder + "reoutput/hist/"
label = "GLISSER"
encounter_file = "encounter.out"


def resonantRatio(a_pl, a_pa):
    klist, jlist, a_rlist = opk.genOuterRes(
        a_pl, 49, 600, high1=5, high2=90, order_lim=90)
    # print(np.sort(a_rlist))

    idx = np.argmin(np.abs(a_pa-a_rlist))

    return klist[idx], jlist[idx], a_rlist[idx]


def plotHistRes(target_pa, target_pl, outputEnc):

    df_pl = pd.read_csv(os.path.join(hist_folder, "pl_{0:d}.txt".format(target_pl)),
                        names=['time', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ')
    # df_pa = pd.read_csv(os.path.join(
    #     histRESO_folder, "pa_{0:d}.csv".format(target_pa)), delimiter=',')
    df_pa = pd.read_csv(os.path.join(
        hist_folder, "pa_{0:d}.txt".format(target_pa)), names=['time', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ')
    df_pl['time'] /= 365*1e6
    df_pl['q'] = df_pl['a'] * (1 - df_pl['e'])
    df_pl['inc'] = np.rad2deg(df_pl['inc'])
    df_pl['Omega'] = np.rad2deg(df_pl['Omega'])
    df_pl['omega'] = np.rad2deg(df_pl['omega'])
    df_pl['M'] = np.rad2deg(df_pl['M'])

    t0, t1 = df_pl['time'].iloc[0], 100

    df_pa['time'] /= 365*1e6
    df_pa['q'] = df_pa['a'] * (1 - df_pa['e'])
    df_pa['inc'] = np.rad2deg(df_pa['inc'])
    df_pa['Omega'] = np.rad2deg(df_pa['Omega'])
    df_pa['omega'] = np.rad2deg(df_pa['omega'])
    df_pa['M'] = np.rad2deg(df_pa['M'])

    df_pa['color'] = opk.BLUE
    df_pa['RESO'] = -2
    df_pa['RESO'] = df_pa['RESO'].shift(periods=182, fill_value=-2)

    df_pa.loc[df_pa['RESO'] > 0, 'color'] = opk.DARK_RED
    df_pa.loc[df_pa['RESO'] < 0, 'color'] = 'gray'
    df_pa.loc[df_pa['RESO'] < -1, 'color'] = 'gray'

    try:
        counts = np.bincount(df_pa['RESO'][df_pa['RESO'] > 0])
        code = np.argmax(counts)
    except:
        code = 0

    # print(df_pl)
    # print(df_pa)

    enc_txt = "reoutput/enc/pa_{0:d}.txt".format(target_pa)
    if outputEnc:
        if (not os.path.exists(folder+enc_txt)):
            print("Retriving close encounter data ...")
            command = "cat {2}{0} | grep ' {1} w/' > {2}reoutput/enc/pa_{1}.txt".format(
                encounter_file, target_pa, folder)
            os.system(command)

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
        time_diff = np.diff(enc_t)
        Period = 5196/1e6
        print(time_diff)
        print(time_diff > Period)
        idxN = (enc_pl == 4)
        idxR = (enc_pl == 5)

    n_pl = df_pl['time'].size
    n_pa = df_pa['time'].size
    print(n_pl, n_pa)

    fig, axes = plt.subplots(4, figsize=(
        8.5, 8), gridspec_kw={'height_ratios': [3, 1, 2, 1]}, sharex=True)
    ax1, ax3, ax2, ax4 = axes[0], axes[1], axes[2], axes[3]
    par_size = 1
    alpha = 1

    ax1.scatter(df_pa['time'], df_pa['a'], alpha=0.5,
                s=par_size*3.5, facecolors=df_pa['color'], edgecolors='none')

    j, k = 0, 0
    if code > 0:
        k = code // 1000
        j = code-k*1000
        loc = opk.resonanceLocationBYA(np.mean(df_pl['a']), k, j)
        ax1.plot([t0, t1], [loc, loc],
                 color='gray', linestyle='dashed', lw=1)

    axE = ax1.twinx()
    if outputEnc and data.shape[0] != 0:
        axE.scatter(enc_t[idxR], enc_dist[idxR], marker='x', s=80, color=opk.DARK_RED, alpha=np.clip(
            enc_dist[idxR]**(-2), a_min=0, a_max=1), lw=2.5)
        axE.scatter(enc_t[idxN], enc_dist[idxN],
                    marker='x', s=50, color='gray')

    ax2.scatter(df_pa['time'], df_pa['q'], alpha=alpha, s=par_size)
    ax5 = ax2.twinx()
    ax5.scatter(df_pa['time'], df_pa['inc'],
                alpha=0.2, s=par_size, color=opk.DARK_RED)

    j = 4
    k = 1
    if code > -3:
        phi = opk.resonantAngleOuter(
            df_pl['Omega'][0:n_pa], df_pl['omega'][0:n_pa], df_pl['M'][0:n_pa], df_pa['Omega'], df_pa['omega'], df_pa['M'], j, k)
        ax3.scatter(df_pa['time'], opk.wrapTo360(phi),
                    alpha=0.5, s=par_size*3.5, color=df_pa['color'], facecolors=df_pa['color'], edgecolors='none')
        loc = opk.resonanceLocationBYA(np.mean(df_pl['a']), k, j)
        ax1.plot([t0, t1], [loc, loc],
                 color='gray', linestyle='dashed', lw=1.5)
        # label = "{0}:{1} Resonant".format(j, k)
        # textcolor = "orange"
        # subsubfolder = "{0}_{1}".format(j, k)
    # if code == 0:
    #     label = "Detached"
    #     textcolor = 'C0'
    #     subsubfolder = "{0}".format(df_pa['reso'])
    # if code < 0:
    #     label = "Scattering"
    #     textcolor = 'red'
    #     subsubfolder = "{0}".format(df_pa['reso'])
    # ax3.text(df_pa['time'][-1]/2, 180, label, color=textcolor,
    #          ha='center', va='center', fontsize=75, weight='bold')

    ax4.scatter(df_pa['time'], opk.wrapTo360(df_pa['omega']),
                alpha=alpha, s=par_size)
    axO = ax4.twinx()
    axO.scatter(df_pa['time'], opk.wrapTo360(df_pa['Omega']),
                alpha=0.2, s=par_size, color=opk.DARK_RED)

    # ax5.scatter(df_pa['time'], opk.wrapTo360(df_pa['omega']+df_pa['Omega']),
    #             alpha=alpha, s=par_size)

    # ax1.text(50, 80, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=36)
    # ax1.text(51, 52.5, "Non-physical", fontsize = 16, color='black',alpha=0.5,rotation=24)

    aN = 30

    # ax2.plot([df_pl['time'][0], df_pl['time'][-1]], [aN, aN], color='red', linestyle='dashed',lw=1)
    ax2.plot([t0, t1], [38, 38], color='gray',
             alpha=1, linestyle='dashed', lw=1.5)
    ax4.plot([t0, t1], [90, 90],
             color='gray', linestyle='dashed', lw=1.5)
    ax4.plot([t0, t1], [270, 270],
             color='gray', linestyle='dashed', lw=1.5)

    ax1.set_title(
        "Particle: {idx:d} - {j:d}:{k:d} Resonance".format(idx=target_pa, j=j, k=k), fontsize=20)
    ax1.set_ylabel("a (au)")
    axE.set_ylabel("Encounter distance (au)", color=opk.DARK_RED)
    axE.set_ylim([0.1, 4])
    axE.set_yscale('log')
    ax2.set_ylabel("q (au)", color=opk.BLUE)
    ax5.set_ylabel("I (deg)", color=opk.DARK_RED)
    ax3.set_ylabel(r"$\varphi_{{{0:d}:{1:d}}}$ (deg)".format(j, k))
    ax3.set_ylim([0, 360])
    ax4.set_ylabel(r"$\omega$ (deg)", color=opk.BLUE)
    axO.set_ylabel(r"$\Omega$ (deg)", color=opk.DARK_RED)
    ax4.set_ylim([0, 360])
    axO.set_ylim([0, 360])

    yvals = ax2.get_yticks()[2:]
    yvals = np.sort(np.append(yvals, 38))
    print(yvals)
    ax2.set_yticks(yvals)
    ax2.set_yticklabels(["{0:d}".format(int(y)) for y in yvals])
    ax2.set_ylim([30, 55])
    # ax5.set_ylabel(r"$\varpi$ (deg)")
    # ax5.set_ylim([0, 360])
    axes[-1].set_xlabel("Time (Myr)", fontsize=18)
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.set_xlim([t0, t1])

    for ax in [ax4, axO, ax3]:
        ax.set_yticks([90, 270])

    axE.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axE.set_yticks([0.1, 0.3, 1, 3])
    ax1.set_ylim(75.1, 76.5)

    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.align_ylabels([ax1, ax2, ax3, ax4])
    fig.align_ylabels([axE, axO, ax5])

    for ax in [axE, axO, ax5]:
        ax.tick_params(axis='y', colors=opk.DARK_RED)

    for ax in [ax2, ax4]:
        ax.tick_params(axis='y', colors=opk.BLUE)

    plt.savefig(output_folder +
                "/particle_{idx:d}_RESO.jpg".format(idx=target_pa), dpi=200)
    print("Saved! particle: {idx:04d}\n".format(idx=target_pa))

    plt.close()

    return None


# target_pl = 4
# output = open(folder + "classes_temp.txt", "w")

for target_pa in [57254]:

    plotHistRes(target_pa, 4, True)

# output.close()
