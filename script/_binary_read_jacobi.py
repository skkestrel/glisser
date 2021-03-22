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




pos_pl, vel_pl, pos_pa, vel_pa = np.empty((0, 3)), np.empty((0, 3)), np.empty((0, 3)), np.empty((0, 3))
fig, ax1 = plt.subplots(figsize=(10,8))

np.set_printoptions(precision=16)
idx=1
for filename,label in zip([filename1],["GLISSER"]):
    t = []
    with open(folder + filename, 'rb') as f:
        read = f.read(24)
        while len(read) == 24:
            time, solar_mass, npl = struct.unpack('=2dQ', read)
            t.append(time)
            for i in range(npl):
                pid, pl_mass, x, y, z, vx, vy, vz = struct.unpack('=I7d', f.read(60)) 
                print(time, x, y, z, vx, vy, vz)
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
