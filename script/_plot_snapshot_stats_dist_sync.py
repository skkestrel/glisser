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

folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Full-Scattering/"
output_folder = folder + "pics/stats/"

for t in [0]:
    a = 50
    q = 32
    rlist = []
    for M in np.linspace(0,360,100000):
        # M = random.uniform(0, 360)
        e = 1-q/a
        f = opk.M2F(np.deg2rad(M), e)
        r = a*(1-e**2)/(1+e*np.cos(f))
        rlist.append(r)
    print(len(rlist))

    x1, x2 = 30, 100
    df = pd.DataFrame(rlist, columns=['r'])
    df = df[df['r'].between(x1, x2)].copy()

    co_list = [-2.5, -3.5]

    xx = np.linspace(x1, x2, 1000)
    fig, ax1 = plt.subplots(1, figsize=(9, 6))
    nbins = 200
    ax1.hist(df['r'], nbins, density=True,color=opk.BLUE, histtype='step')
    for co, color in zip(co_list, [opk.GREEN, opk.RED]):
        scale = -1/(1+co) * x1**(1+co) + 1/(1+co) * x2**(1+co)

        yy = xx**(co)/scale


        # ax1.plot(xx, yy, color=color, ls='dashed',
        #          label="co = {0:.1f}".format(co))

    # a_N = 30
    # klist, jlist, a_rlist = opk.genOuterRes(
    #     a_N, 50, 125, high1=2, high2=50, order_lim=11)
    # for k, j, a_res in zip(klist, jlist, a_rlist):
    #     label = "{0}/{1}".format(j, k)
    #     ax1.text(a_res-0.3, 0.065, label, fontsize=10, rotation=70)
    #     # print(a)
    #     ax1.plot([a_res, a_res], [0, 300], ls='dashed',
    #              color='gray', zorder=1, alpha=0.6)

    ax1.set_xlabel("a (au)")
    ax1.set_ylabel("Count")
    ax1.set_xlim(x1-5, x2+5)
    ax1.set_ylim(0, 0.4)
    ax1.legend()
    plt.title("a={0:d}, q={1:d}".format(a, q))
    plt.tight_layout()

    plt.savefig(output_folder +
                "sythetic_a={0:d}, q={1:d}.jpg".format(a, q), dpi=200)
    print("Saved! frame")
    plt.close()
