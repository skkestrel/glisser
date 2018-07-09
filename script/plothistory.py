import sys
import math
import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import struct
import matplotlib.style as style

style.use('ggplot')

times = []
planets = [[] for i in range(12)]
particlewatch = [int(sys.argv[i + 3]) for i in range(len(sys.argv) - 3)]
particles = [[] for i in range(len(particlewatch) * 3)]

def bsearch(f, npa, partnum):
    base = f.tell()
    left = 0
    right = npa-1

    while left <= right:
        mid = left + (right - left) // 2
        f.seek(base + mid * 16, 0)
        Id, a, e, i = struct.unpack('=I3f', f.read(16))

        if Id == partnum:
            f.seek(base, 0)
            return Id, a, e, i
        elif Id > partnum:
            right = mid - 1
        else:
            left = mid + 1

    f.seek(base, 0)
    return None

take_every = 1
if len(sys.argv) > 2:
    take_every = int(sys.argv[2])

with open(sys.argv[1], 'rb') as f:
    read = f.read(16)
    counter = 0
    while len(read) == 16:
        time, npl = struct.unpack('=dQ', read)
        if (counter % take_every) == 0:
            times.append(time)

        for i in range(npl):
            a, e, I = struct.unpack('=I3f', f.read(16))[1:]

            if (counter % take_every) == 0:
                planets[3*i].append(a)
                planets[3*i+1].append(e)
                planets[3*i+2].append(I)
        npa, = struct.unpack('=Q', f.read(8))

        if (counter % take_every) == 0:
            for index, partnum in enumerate(particlewatch):
                part = bsearch(f, npa, partnum)
                if part:
                    particles[3*index].append(part[1])
                    particles[3*index+1].append(part[2])
                    particles[3*index+2].append(part[3])
                else:
                    particles[3*index].append(math.nan)
                    particles[3*index+1].append(math.nan)
                    particles[3*index+2].append(math.nan)

        f.seek(npa * 16, 1)
        read = f.read(16)
        counter = counter + 1;

times = np.array(times)

planets = np.array(planets)
particles = np.array(particles)
cc = plt.rcParams['axes.prop_cycle'].by_key()['color']


plt.figure()

if times[-1] > 365e6:
    stimes = times / 365e6
    timelabel = "Time (Myr)"
elif times[-1] > 365e3:
    stimes = times / 365e3
    timelabel = "Time (kyr)"
else:
    stimes = times
    timelabel = "Time (yr)"

plt.xlabel(timelabel)

for i in range(4):
    c = cc[(i) % len(cc)]
    plt.plot(stimes, planets[3*i], c=c)
    plt.plot(stimes, planets[3*i] * (1. - planets[3*i+1] ), c=c)
    plt.plot(stimes, planets[3*i] * (1. + planets[3*i+1] ), c=c)
    plt.plot([], [], c=cc[i%len(cc)], label="Planet {0}".format(i+1))
for i in range(len(particlewatch)):
    c = cc[(i+4) % len(cc)]
    nnull = particles[2*i] != math.nan

    plt.plot(stimes[nnull], particles[3*i][nnull], c=c)
    plt.plot(stimes[nnull], (particles[3*i] * ( 1. - particles[3*i+1]) )[nnull], c=c)
    plt.plot(stimes[nnull], (particles[3*i] * ( 1. + particles[3*i+1]) )[nnull], c=c)
    plt.plot([], [], c=c, label="Particle {0}".format(particlewatch[i]))

plt.title(sys.argv[1])
plt.legend()

fig, axes = plt.subplots(3, 1, sharex=True)
dt = times[1] - times[0]
taxis = np.fft.fftfreq(len(times), dt / 365) * 360 * 3600
half = len(times) // 2

for i in range(4):
    c = cc[(i) % len(cc)]
    axes[0].plot(taxis[:half], np.abs(np.fft.fft(planets[3*i])[:half]), c=c)
    axes[1].plot(taxis[:half], np.abs(np.fft.fft(planets[3*i+1])[:half]), c=c)
    axes[2].plot(taxis[:half], np.abs(np.fft.fft(planets[3*i+2])[:half]), c=c)
    axes[0].plot([], [], c=cc[i%len(cc)], label="Planet {0}".format(i+1))
axes[0].set_title("a")
axes[0].set_yscale("log")
axes[0].set_xscale("log")
axes[1].set_title("e")
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[2].set_title("i")
axes[2].set_xscale("log")
axes[2].set_yscale("log")
axes[2].set_xlabel("arcsec / yr")
axes[0].legend()


fig, axes = plt.subplots(3, 1, sharex=True)
for i in range(4):
    c = cc[(i) % len(cc)]
    axes[0].plot(stimes, planets[3*i]-planets[3*i].mean(), c=c)
    axes[1].plot(stimes, planets[3*i+1], c=c)
    axes[2].plot(stimes, planets[3*i+2], c=c)
    axes[0].plot([], [], c=cc[i%len(cc)], label="Planet {0}".format(i+1))
axes[0].set_title("Δa")
axes[1].set_title("e")
axes[2].set_title("i")
axes[2].set_xlabel(timelabel)
axes[0].legend()
plt.show()
