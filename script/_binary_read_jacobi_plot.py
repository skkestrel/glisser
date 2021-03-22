import numpy as np
import matplotlib.pyplot as plt
import struct

folder = "/home/yhuang/GLISSER/glisser/"
filename1 = "benchmark/circular-neptune-particles-test/tracks/track.0.out"
discard_file = "benchmark/circular-neptune-particles-test/encounter.out"


def jacobi_integral(solar_mass, pl_mass, omega, pos_pl, vel_pl, pos_pa, vel_pa):
    pos_plh, vel_plh, pos_sun, vel_sun = coord_b2h(solar_mass, pl_mass, pos_pl, vel_pl)

    jacobi = 0.5 * (vel_pa[:,0]**2 + vel_pa[:,1]**2 + vel_pa[:,2]**2) 

    rr = np.sqrt((pos_pa[:,0] - pos_sun[:,0])**2 \
                + (pos_pa[:,1] - pos_sun[:,1])**2 \
                + (pos_pa[:,2] - pos_sun[:,2])**2)

    jacobi -= solar_mass/rr

    rr = np.sqrt((pos_pa[:,0] - pos_pl[:,0])**2 \
                + (pos_pa[:,1] - pos_pl[:,1])**2 \
                + (pos_pa[:,2] - pos_pl[:,2])**2)

    jacobi -= pl_mass/rr

    jx = pos_pa[:,1] * vel_pa[:,2] - pos_pa[:,2] * vel_pa[:,1]
    jy = pos_pa[:,2] * vel_pa[:,0] - pos_pa[:,0] * vel_pa[:,2]
    jz = pos_pa[:,0] * vel_pa[:,1] - pos_pa[:,1] * vel_pa[:,0]

    jacobi -= (omega[:,0] * jx + omega[:,1] * jy + omega[:,2] * jz)

    return jacobi

def get_omega(solar_mass, pl_mass, pos_plb, vel_plb):
    pos_pl, vel_pl, pos_sun, vel_sun = coord_b2h(solar_mass, pl_mass, pos_plb, vel_plb)
    gm = solar_mass + pl_mass
    energy = 0.5 * (vel_pl[:,0]**2 + vel_pl[:,1]**2 + vel_pl[:,2]**2)
    energy -= gm/np.sqrt(pos_pl[:,0]**2 + pos_pl[:,1]**2 + pos_pl[:,2]**2)
    aplh = -0.5 * gm/energy
    omega = np.sqrt(gm/(aplh**3))
    omegax = pos_pl[:,1] * vel_pl[:,2] - pos_pl[:,2] * vel_pl[:,1]
    omegay = pos_pl[:,2] * vel_pl[:,0] - pos_pl[:,0] * vel_pl[:,2]
    omegaz = pos_pl[:,0] * vel_pl[:,1] - pos_pl[:,1] * vel_pl[:,0]
    fac = omega / np.sqrt(omegax**2 + omegay**2 + omegaz**2)
    omegax = omegax * fac
    omegay = omegay * fac
    omegaz = omegaz * fac
    return np.hstack((omegax.reshape(omegax.shape[0],-1), \
                      omegay.reshape(omegay.shape[0],-1), \
                      omegaz.reshape(omegaz.shape[0],-1)))


def coord_b2h(solar_mass, pl_mass, pos_plb, vel_plb):
    ratio = pl_mass/(pl_mass + solar_mass)
    pos_pl = pos_plb/(1-ratio)
    vel_pl = vel_plb/(1-ratio)

    center_r = pos_pl - pos_plb
    center_v = vel_pl - vel_plb
    return pos_pl, vel_pl, -center_r, -center_v



files = [filename1]


# fig2, ax6 = plt.subplots(figsize=(9,6))
np.set_printoptions(precision=16)

discard_list = np.loadtxt(discard_file, usecols=(1, 5), delimiter=' ')
discard_idx = [1,2,3]

for discard in discard_list:
    idx = int(discard[0])
    if(idx in discard_idx):
        continue

    discard_idx.append(int(discard[0]))
    pos_pl, vel_pl, pos_pa, vel_pa = np.empty((0, 3)), np.empty((0, 3)), np.empty((0, 3)), np.empty((0, 3))
    fig, ax1 = plt.subplots(figsize=(10,8))
    for filename,label in zip(files,["GLISSER"]):
        t = []
        with open(folder + filename, 'rb') as f:
            read = f.read(24)
            while len(read) == 24:
                time, solar_mass, npl = struct.unpack('=2dQ', read)
                t.append(time)
                for i in range(npl):
                    pid, pl_mass, x, y, z, vx, vy, vz = struct.unpack('=I7d', f.read(60)) 
                    pos_pl = np.vstack((pos_pl, np.array([x, y, z])))
                    vel_pl = np.vstack((vel_pl, np.array([vx, vy, vz])))
                    
                npa, = struct.unpack('=Q', f.read(8))
                # print(npa)
                for i in range(npa):
                    pid, x, y, z, vx, vy, vz = struct.unpack('=I6d', f.read(52))
                    if(pid==idx):
                        pos_pa = np.vstack((pos_pa, np.array([x, y, z])))
                        vel_pa = np.vstack((vel_pa, np.array([vx, vy, vz])))
                read = f.read(24)
        num = len(t)
        
        omega = get_omega(solar_mass, pl_mass, pos_pl, vel_pl)
        Cj = jacobi_integral(solar_mass, pl_mass, omega, pos_pl, vel_pl, pos_pa, vel_pa)


        
        size = 10
        ax1.scatter(t, np.abs((Cj-Cj[0])/Cj[0]), alpha=0.5, label=label, s=size)
        ax1.grid(linestyle='dashed',alpha=0.3)
        # ax1.plot([2.5e6, 2.5e6], [1e-14,1e-2], ls='dashed',color='gray', alpha=0.5,lw=1)
        # ax1.plot([5e6, 5e6], [1e-14,1e-2], ls='dashed',color='gray',alpha=0.5, lw=1)
        # ax1.plot([7.5e6, 7.5e6], [1e-14,1e-2], ls='dashed',color='gray', alpha=0.5,lw=1)
        for discard in discard_list:
            if(abs(idx-discard[0])<1e-5):
                discard_time = discard[1]
                ax1.plot([discard_time , discard_time], [1e-14,1e-5] ,color='red', alpha=0.5,lw=1)

    ax1.set_xlabel("Time (days)")
    ax1.set_xlabel("Particles {id} (stepsize = 300 days)".format(id=idx))   
    ax1.set_yscale("log")
    ax1.set_ylim(1e-14,1e-5)
    # ax2.legend()
    print("Saving: {id}".format(id=idx))
    fig.savefig("jacobi_300/jacobi_erros_{id}.jpg".format(id=idx),dpi=200)
    plt.close()
    # fig2.savefig("tisserand3_particle_{id}_Cj_comp2.png".format(id=idx),dpi=200)