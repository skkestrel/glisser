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
particlewatch = [int(sys.argv[i + 6]) for i in range(len(sys.argv) - 6)]
particles = [[] for i in range(len(particlewatch) * 6)]

def bsearch(f, npa, partnum):
    base = f.tell()
    left = 0
    right = npa-1

    while left <= right:
        mid = left + (right - left) // 2
        f.seek(base + mid * 28, 0)
        Id, a, e, i, O, o, f = struct.unpack('=I6f', f.read(28))

        if Id == partnum:
            f.seek(base, 0)
            return Id, a, e, i, O, o, f
        elif Id > partnum:
            right = mid - 1
        else:
            left = mid + 1

    f.seek(base, 0)
    return None

take_every = 1
if len(sys.argv) > 2:
    take_every = int(sys.argv[2])

planets = None

with open(sys.argv[1], 'rb') as f:
    read = f.read(16)
    counter = 0
    while len(read) == 16:
        time, npl = struct.unpack('=dQ', read)

        if not planets:
            planets = list([[] for i in range(npl * 6)])

        if (counter % take_every) == 0:
            times.append(time)

        for i in range(npl):
            a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))[1:]

            if (counter % take_every) == 0:
                planets[6*i].append(a)
                planets[6*i+1].append(e)
                planets[6*i+2].append(I)
                planets[6*i+3].append(O)
                planets[6*i+4].append(o)
                planets[6*i+5].append(F)
        npa, = struct.unpack('=Q', f.read(8))

        if (counter % take_every) == 0:
            for index, partnum in enumerate(particlewatch):
                part = bsearch(f, npa, partnum)
                if part:
                    particles[6*index].append(part[1])
                    particles[6*index+1].append(part[2])
                    particles[6*index+2].append(part[6])
                    particles[6*index+6].append(part[4])
                    particles[6*index+4].append(part[5])
                    particles[6*index+5].append(part[6])
                else:
                    particles[6*index].append(math.nan)
                    particles[6*index+1].append(math.nan)
                    particles[6*index+2].append(math.nan)
                    particles[6*index+6].append(math.nan)
                    particles[6*index+4].append(math.nan)
                    particles[6*index+5].append(math.nan)

        f.seek(npa * 28, 1)
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

for i in range(npl):
    c = cc[(i) % len(cc)]
    plt.plot(stimes, planets[6*i], c=c)
    plt.plot(stimes, planets[6*i] * (1. - planets[6*i+1] ), c=c)
    plt.plot(stimes, planets[6*i] * (1. + planets[6*i+1] ), c=c)
    plt.plot([], [], c=cc[i%len(cc)], label="Planet {0}".format(i+1))
for i in range(len(particlewatch)):
    c = cc[(i+npl) % len(cc)]
    nnull = particles[6*i] != math.nan

    plt.plot(stimes[nnull], particles[6*i][nnull], c=c)
    plt.plot(stimes[nnull], (particles[6*i] * ( 1. - particles[6*i+1]) )[nnull], c=c)
    plt.plot(stimes[nnull], (particles[6*i] * ( 1. + particles[6*i+1]) )[nnull], c=c)
    plt.plot([], [], c=c, label="Particle {0}".format(particlewatch[i]))

plt.title(sys.argv[1])
plt.legend()




fig, axes = plt.subplots(3, 1, sharex=True)
dt = times[1] - times[0]
taxis = np.fft.fftfreq(len(times), dt / 365) * 360 * 3600
half = len(times) // 2

for i in range(npl):
    c = cc[(i) % len(cc)]
    axes[0].plot(taxis[1:half], np.abs(np.fft.fft(planets[6*i] - planets[6*i].mean())[1:half]), c=c)
    axes[1].plot(taxis[1:half], np.abs(np.fft.fft(planets[6*i+1] - planets[6*i+1].mean())[1:half]), c=c)
    axes[2].plot(taxis[1:half], np.abs(np.fft.fft(planets[6*i+2] - planets[6*i+2].mean())[1:half]), c=c)
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




fig, axes = plt.subplots(2, 1, sharex=True)
for i in range(npl):
    c = cc[(i) % len(cc)]
    axes[0].plot(taxis[1:half], np.abs(np.fft.fft(planets[6*i+2] * np.sin(planets[6*i+3]))[1:half]), c=c)
    axes[1].plot(taxis[1:half], np.abs(np.fft.fft(planets[6*i+1] * np.sin(planets[6*i+4] + planets[6*i+3]))[1:half]), c=c)
    axes[0].plot([], [], c=cc[i%len(cc)], label="Planet {0}".format(i+1))
axes[0].set_title("p")
axes[0].set_yscale("log")
axes[0].set_xscale("log")
axes[1].set_title("h")
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_xlabel("arcsec / yr")
axes[0].legend()




fig, axes = plt.subplots(3, 1, sharex=True)
for i in range(npl):
    c = cc[(i) % len(cc)]
    axes[0].plot(stimes, planets[6*i]-planets[6*i].mean(), c=c)
    axes[1].plot(stimes, planets[6*i+1], c=c)
    axes[2].plot(stimes, planets[6*i+2], c=c)
    axes[0].plot([], [], c=cc[i%len(cc)], label="planet {0}".format(i+1))
axes[0].set_title("δa")
axes[1].set_title("e")
axes[2].set_title("i")
axes[2].set_xlabel(timelabel)
axes[0].legend()



def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
fig, axes = plt.subplots(1, 1, sharex=True)
for i in range(npl):
    c = cc[(i) % len(cc)]
    axes.plot(stimes[9:], moving_average(planets[6*i+1], 10), c=c)
    axes.plot([], [], c=cc[i%len(cc)], label="planet {0}".format(i+1))
axes.set_title("e")
axes.set_xlabel(timelabel)
axes.legend()





fig, axes = plt.subplots(3, 1, sharex=True)
for i in range(npl):
    F = planets[6*i+5].copy()

    def normalize(x):
        X = x.copy()

        diff = 0;
        for j in range(1, len(X)):
            X[j] += diff
            if X[j - 1] - X[j] > math.pi:
                X[j] += 2 * math.pi;
                diff += 2 * math.pi;
            elif X[j - 1] - X[j] < -math.pi:
                X[j] -= 2 * math.pi;
                diff -= 2 * math.pi;
        import scipy
        return np.gradient(X, dt / 365.25)

    c = cc[(i) % len(cc)]
    axes[0].plot(stimes, normalize(planets[6*i+3]), c=c)
    axes[1].plot(stimes, normalize(planets[6*i+3]+planets[6*i+4]), c=c)
    axes[2].plot(stimes, normalize(planets[6*i+5]), c=c)
    axes[0].plot([], [], c=cc[i%len(cc)], label="planet {0}".format(i+1))

axes[0].set_title("dΩ/dt (rad/yr)")
axes[1].set_title("d(Ω+ω)/dt")
axes[2].set_title("df/dt")
axes[2].set_xlabel(timelabel)
axes[0].legend()


plt.show()
