import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import pandas as pd
import matplotlib
import orbitplotkit as opk
import random

font = {'weight': 'bold',
        'size': 14}
matplotlib.rc('font', **font)
light_blue = '#5ed2f1'
dark_blue = '#0084c6'
yellow = '#ffb545'
red = '#ff2521'


# folder = "/home/yhuang/GLISSER/glisser/fast/cold_output/out-ivp-ifree-low/reoutput/"
folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-N-scattering-3/"
output_folder = folder + "pics/power_law/"
label = "GLISSER"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# time = [0, 182000000000]
# time = [6275000000]
# time1 = np.arange(4000000000, 73100000001, 1000000000)
time = np.arange(0, 182400000001, 1000000000)
# time = np.concatenate([time2, time1])
# print(time.shape[0])
idx = 0
# time = [182400000000]
# time = [0]
par_size = 15

pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(0)
df_pa0 = opk.snapshot2df(folder+pa_txt)

for t in time:
    pl_txt = "reoutput/snapshots/planets_{0:d}.txt".format(t)
    df_pl = opk.snapshot2df(folder+pl_txt)
    aN = df_pl['a'].iloc[-1]


    print(aN)
    pa_txt = "reoutput/snapshots/particles_{0:d}.txt".format(t)
    df_pa = opk.snapshot2df(folder+pa_txt,apl=aN)
    # print(df_pa)

    x1, x2 = df_pa['a'].min(), 2000
    x1 = 40
    # df_pa = df_pa[df_pa['id'] <= 200000]
    df_pa = df_pa[df_pa['a'].between(x1, x2)]
    # df_pa = df_pa[df_pa['q'] < 40]

    # df_pa['vinf'] = opk.vinf(df_pa['a'], df_pa['e'],
    #                          np.deg2rad(df_pa['inc']), aN)
    # df_pa['vpa'] = opk.vpa(df_pa['a'], aN)
    # df_pa['alpha'] = opk.alpha(df_pa['vinf'], df_pa['vpa'])
    print(df_pa)
    # df_pa['alpha']
    # rlist = []
    # for index, row in df_pa.iterrows():
    #     for i in np.arange(10):
    #         M = random.uniform(0, 360)
    #         a = row['a']
    #         e = row['e']
    #         f = opk.M2F(np.deg2rad(M), e)
    #         r = a*(1-e**2)/(1+e*np.cos(f))
    #         rlist.append(r)
    # print(len(rlist))

    # df = pd.DataFrame(rlist, columns=['r'])
    # df = df[df['a'].between(x1, x2)].copy()
    
    co_list = [-1.5, -2.5, -3.5]
    xx = np.linspace(x1, x2, 2000)
    nbins = 100
    logbins = np.logspace(np.log10(x1),np.log10(x2),nbins)
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    ax1, ax2, ax3, ax4 = axes[0][0], axes[1][0], axes[0][1], axes[1][1]
    ax1.hist(df_pa['a'], bins=logbins, density=True,
                 color=opk.BLUE, histtype='step', log=True)

    for co, color in zip(co_list, [opk.RED, opk.GREEN, opk.LIGHT_BLUE]):
        scale = -1/(1+co) * x1**(1+co) + 1/(1+co) * x2**(1+co)
        yy = xx**(co)/scale
        ax1.plot(xx, yy, color=color, ls='dashed',
                 label="co = {0:.1f}".format(co))

    ax1.set_xlabel("a (au)")
    ax1.set_ylabel("Count")
    ax1.set_xlim(x1-5, x2+5)
    ax1.set_ylim(1e-6, 1e-1)
    ax1.set(xscale='log', yscale='log')
    ax1.legend(loc='upper right',fontsize=12)

    ax2.hist(df_pa['q'], nbins, density=True,
             color=opk.BLUE, histtype='step')
    # vmax = np.sqrt(2)-1
    # ax2.plot([vmax, vmax], [0, 3], ls='dashed', c=opk.DARK_RED)
    ax2.set_xlim(10, 40)
    ax2.set_ylim(0, 1)
    ax2.set_xlabel("q (au)")

    ax3.hist(df_pa['inc'], nbins, density=True,
             color=opk.BLUE, histtype='step')
    ax3.set_xlim(0, 30)
    ax3.set_ylim(0, 1)
    ax3.set_xlabel("i (deg)")

    ax4.hist(df_pa['T'], nbins, density=True,
             color=opk.BLUE, histtype='step')
    ax4.set_xlim(2.95, 3.05)
    ax4.set_ylim(0, 1000)
    ax4.set_xlabel("T_n")
    # ax3.set_xlabel("alpha")
    # ax3.set_xlim(0, np.pi)
    # xx = np.linspace(0, np.pi/2,1000)
    # yy = np.sin(xx)
    # ax3.plot(xx, yy, ls='dashed',c=opk.DARK_RED)


    plt.title("Time: {time:6.3f} Myr".format(time=t/365.25/1e6))
    plt.tight_layout()

    plt.savefig(output_folder +
                "power_law_{idx:04d}.jpg".format(idx=idx), dpi=200)
    print("Saved! frame: {idx:04d}".format(idx=idx))
    plt.close()
    idx += 1
