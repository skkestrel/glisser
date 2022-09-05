import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import os
import matplotlib
import rebound
import orbitplotkit as opk


font = {'weight': 'bold',
        'size': 16}
matplotlib.rc('font', **font)


folder = "/home/yhuang/GLISSER/glisser/fast/rogue_output/out-JSUNT-Giant-Synthetic/"
output_folder = folder + "pics/hist/3D_animation/"
enc_folder = folder + "reoutput/enc/"
histRESO_folder = folder + "reoutput/hist_RESO/"
hist_folder = folder + "reoutput/hist/"
label = "GLISSER"
encounter_file = "encounter.out"


def plotOrbit(target_pa, target_pl, time, showLocations):
    num = time.shape[0]
    print(num)
    idx = 1
    for t in time:
        df_pl = pd.read_csv(os.path.join(hist_folder, "pl_{0:d}.txt".format(target_pl)),
                            names=['time', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ')
        # df_pa = pd.read_csv(os.path.join(
        #     histRESO_folder, "pa_{0:d}.csv".format(target_pa)), delimiter=',')
        df_pa = pd.read_csv(os.path.join(
            hist_folder, "pa_{0:d}.txt".format(target_pa)), names=['time', 'a', 'e', 'inc', 'Omega', 'omega', 'M'], delimiter=' ')
        # df_pl['time'] /= 365*1e6
        pl_orbit = df_pl[df_pl['time'] == t]
        pa_orbit = df_pa[df_pa['time'] == t]
        # print(pl_orbit['a'].values[0])

        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=0, a=pl_orbit['a'].values[0], e=pl_orbit['e'].values[0], inc=pl_orbit['inc'].values[0],
                Omega=pl_orbit['Omega'].values[0], omega=pl_orbit['omega'].values[0], M=pl_orbit['M'].values[0])
        sim.add(m=0, a=pa_orbit['a'].values[0], e=pa_orbit['e'].values[0], inc=pa_orbit['inc'].values[0],
                Omega=pa_orbit['Omega'].values[0], omega=pa_orbit['omega'].values[0], M=pa_orbit['M'].values[0])

        npts = 256
        pl = sim.particles[1]
        pts = np.array(pl.sample_orbit(
            Npts=npts, primary=sim.particles[0]))/100
        pl_x, pl_y, pl_z = [pts[:, i] for i in range(3)]
        node_crossings = np.where(np.diff(np.sign(pl_z)))[0]
        if node_crossings.shape[0] < 2:
            node_crossings = np.append(node_crossings, -1)
        df_pl = pd.DataFrame({'x': pl_x, 'y': pl_y, 'z': pl_z})

        if df_pl['z'].iloc[0] > 0:
            df_pl1 = pd.concat(
                (df_pl[node_crossings[1]+1:], df_pl[:node_crossings[0]+1]))
            df_pl2 = df_pl[node_crossings[0]+1:node_crossings[1]+1]
        else:
            df_pl2 = pd.concat(
                (df_pl[node_crossings[1]+1:], df_pl[:node_crossings[0]+1]))
            df_pl1 = df_pl[node_crossings[0]+1:node_crossings[1]+1]
        df_pl_node = df_pl.iloc[node_crossings]
        # print(df_pl1, df_pl2)

        pa = sim.particles[2]
        pts = np.array(pa.sample_orbit(
            Npts=npts, primary=sim.particles[0]))/100
        pa_x, pa_y, pa_z = [pts[:, i] for i in range(3)]
        node_crossings = np.where(np.diff(np.sign(pa_z)))[0]
        if node_crossings.shape[0] < 2:
            node_crossings = np.append(node_crossings, -1)
        # print(node_crossings, pa_z[node_crossings])

        df_pa = pd.DataFrame({'x': pa_x, 'y': pa_y, 'z': pa_z})
        if df_pa['z'].iloc[0] > 0:
            df_pa1 = pd.concat(
                (df_pa[node_crossings[1]+1:], df_pa[:node_crossings[0]+1]))
            df_pa2 = df_pa[node_crossings[0]+1:node_crossings[1]+1]
        else:
            df_pa2 = pd.concat(
                (df_pa[node_crossings[1]+1:], df_pa[:node_crossings[0]+1]))
            df_pa1 = df_pa[node_crossings[0]+1:node_crossings[1]+1]
        df_pa_node = df_pa.iloc[node_crossings]
        # print(df_pa_node)
        # print(df_pl)

        fig = plt.figure(figsize=(8, 8))

        ax = plt.axes(projection='3d', computed_zorder=False)
        angle = (idx-1)/num * 360
        ax.view_init(30, angle)
        print(angle)

        xx, yy = np.meshgrid(np.linspace(-1,1), np.linspace(-1,1))
        zz = yy*0
        ax.plot_surface(xx, yy, zz, alpha=0.5, color=opk.GRAY)

        ax.plot3D(df_pl1['x'], df_pl1['y'],
                  df_pl1['z'], c=opk.RED, lw=2, zorder=2)
        ax.plot3D(df_pl2['x'], df_pl2['y'], df_pl2['z'],
                  c=opk.RED, alpha=0.25, lw=2, zorder=2)
        ax.scatter(df_pl_node['x'], df_pl_node['y'], df_pl_node['z'], alpha=1,
                   marker='.', s=100, facecolors=opk.RED, edgecolor='none', zorder=4)

        ax.plot3D(df_pa1['x'], df_pa1['y'],
                  df_pa1['z'], c=opk.LIGHT_BLUE, lw=2)
        ax.plot3D(df_pa2['x'], df_pa2['y'], df_pa2['z'],
                  c=opk.LIGHT_BLUE, alpha=0.25, lw=2)
        ax.scatter(df_pa_node['x'], df_pa_node['y'], df_pa_node['z'], alpha=1,
                   marker='.', s=100, facecolors=opk.LIGHT_BLUE, edgecolor='none', zorder=4)

        ax.grid(ls='dashed', alpha=0.1, lw=0.5)
        ax.set(xlim=[-1, 1], ylim=[-1, 1], zlim=[-0.5, 0.5],
               xticks=[-1, 0, 1], yticks=[-1, 0, 1], zticks=[-0.5, 0, 0.5])

        ax.set_title("Time: {0:.2f} Myr".format(t/365.25/1e6), fontsize=20)
        # fig[0].tight_layout()
        plt.savefig(output_folder +
                    "/particle_{pa:d}_Orbit_{idx:04d}.jpg".format(idx=idx, pa=target_pa), dpi=200)
        print("Saved! particle: {pa:04d}\n".format(pa=target_pa))

        plt.close()
        idx += 1

    return None


# target_pl = 4
# output = open(folder + "classes_temp.txt", "w")

for target_pa in [57254]:
    time = np.arange(0, 36520000001, 50000000)
    plotOrbit(target_pa, 5, time, False)

# output.close()
