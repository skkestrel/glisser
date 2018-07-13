"""plot_history.py

Usage:
    plot_history.py [options] <path>

Options:
    -h, --help                     Show this screen.
    -w <list>, --watch <list>      Plot the comma-separated list of particles, or "all"
    -P, --no-planets               Don't plot planets.
    -s <n>, --skip <n>             Take every n time steps [default: 1]
    -t <t>, --tmax <t>             Take only up to given time
    --plot-angles                  ..
    --plot-aei                     ..
    --plot-ftrect                  ..
    --plot-qQ                      ..
    --plot-qQ-mmr                  ..
    --plot-e-smooth                ..
"""

import sys
import os
import math
import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import struct
import matplotlib.style as style
import docopt

args = docopt.docopt(__doc__)

style.use('ggplot')

times = []

if args["--watch"]:
    if args["--watch"] == "all":
        particlewatch = []
        particlewatch_rev = {}
        take_all_particles = True
    else:
        take_all_particles = False
        particlewatch = [int(x) for x in args["--watch"].split(',')]
else:
    take_all_particles = False
    particlewatch = []

particles = [[] for i in range(len(particlewatch) * 6)]

def bsearch(f, npa, partnum):
    base = f.tell()
    left = 0
    right = npa-1

    while left <= right:
        mid = left + (right - left) // 2
        f.seek(base + mid * 28, 0)
        Id, a, e, i, O, o, F = struct.unpack('=I6f', f.read(28))

        if Id == partnum:
            f.seek(base, 0)
            return Id, a, e, i, O, o, F
        elif Id > partnum:
            right = mid - 1
        else:
            left = mid + 1

    f.seek(base, 0)
    return None

take_every = int(args["--skip"])

planets = None
skip_planets = bool(args["--no-planets"])
if skip_planets:
    planets = []

filenum = 0
counter = 0
while True:
    try:
        if os.path.isfile(args["<path>"]):
            filename = args["<path>"]
        elif os.path.isdir(args["<path>"]):
            filename = os.path.join(args["<path>"], "track.{0}.out".format(filenum))
        else:
            print("File does not exist")
            sys.exit(-1)

        with open(filename, 'rb') as f:
            filenum += 1
            read = f.read(16)

            while len(read) == 16:
                time, npl = struct.unpack('=dQ', read)

                if not planets and not skip_planets:
                    planets = list([[] for i in range(npl * 6)])

                if args["--tmax"] and time > float(args["--tmax"]):
                    raise IOError() # dirty, but just break out

                if (counter % take_every) == 0:
                    times.append(time)

                if (counter % take_every) == 0 and not skip_planets:
                    for i in range(npl):
                        a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))[1:]

                        planets[6*i].append(a)
                        planets[6*i+1].append(e)
                        planets[6*i+2].append(I)
                        planets[6*i+3].append(O)
                        planets[6*i+4].append(o)
                        planets[6*i+5].append(F)
                else:
                    f.seek(npl * 28, 1)

                npa, = struct.unpack('=Q', f.read(8))

                if (counter % take_every) == 0:
                    if take_all_particles:
                        if len(particlewatch) == 0:
                            firstrun = True
                        else:
                            firstrun = False

                        foundparticles = set()

                        for i in range(npa):
                            pid, a, e, I, O, o, F = struct.unpack('=I6f', f.read(28))

                            if firstrun:
                                particlewatch_rev[pid] = len(particlewatch)
                                particlewatch.append(pid)
                                for j in range(6):
                                    particles.append([])

                            foundparticles.add(pid)
                            ind = particlewatch_rev[pid]


                            particles[6*ind].append(a)
                            particles[6*ind+1].append(e)
                            particles[6*ind+2].append(I)
                            particles[6*ind+3].append(O)
                            particles[6*ind+4].append(o)
                            particles[6*ind+5].append(F)

                        for i in set(particlewatch) - foundparticles:
                            ind = particlewatch_rev[i]

                            particles[6*ind].append(math.nan)
                            particles[6*ind+1].append(math.nan)
                            particles[6*ind+2].append(math.nan)
                            particles[6*ind+3].append(math.nan)
                            particles[6*ind+4].append(math.nan)
                            particles[6*ind+5].append(math.nan)

                    else:
                        for index, partnum in enumerate(particlewatch):
                            part = bsearch(f, npa, partnum)
                            if part:
                                particles[6*index].append(part[1])
                                particles[6*index+1].append(part[2])
                                particles[6*index+2].append(part[3])
                                particles[6*index+3].append(part[4])
                                particles[6*index+4].append(part[5])
                                particles[6*index+5].append(part[6])
                            else:
                                particles[6*index].append(math.nan)
                                particles[6*index+1].append(math.nan)
                                particles[6*index+2].append(math.nan)
                                particles[6*index+3].append(math.nan)
                                particles[6*index+4].append(math.nan)
                                particles[6*index+5].append(math.nan)
                        f.seek(npa * 28, 1)
                else:
                    f.seek(npa * 28, 1)

                read = f.read(16)
                counter = counter + 1;
        if os.path.isfile(args["<path>"]):
            break
    except IOError:
        break



times = np.array(times)
planets = np.array(planets)
particles = np.array(particles)

cc = plt.rcParams['axes.prop_cycle'].by_key()['color']

if times[-1] > 365e6:
    stimes = times / 365e6
    timelabel = "Time (Myr)"
elif times[-1] > 365e3:
    stimes = times / 365e3
    timelabel = "Time (kyr)"
else:
    stimes = times
    timelabel = "Time (yr)"

if skip_planets:
    npl = 0

def do_for(callback):
    for i in range(npl):
        c = cc[i % len(cc)]
        callback(planets[6*i:6*i+6], c, i+1, True)
    for i in range(len(particlewatch)):
        c = cc[(i+npl) % len(cc)]
        callback(particles[6*i:6*i+6], c, particlewatch[i], False)

if args["--plot-qQ"]:
    plt.figure()
    def plot_qQ(data, c, ind, is_planet):
        plt.plot(stimes, data[0], c=c)
        plt.plot(stimes, data[0] * (1. - data[1]), c=c)
        plt.plot(stimes, data[0] * (1. + data[1]), c=c)
        if is_planet:
            if args["--plot-qQ-mmr"] and ind != npl:
                def mmr(a1, deg):
                    return (deg * (a1 ** (3/2))) ** (2/3)
                
                meana = data[0].mean()
                plt.axhline(y=meana, linewidth=2, color=c)

                for i in range(1, 5):
                    res = mmr(meana, (i+1)/i)
                    plt.axhline(y=res, linewidth=1, color=c, zorder=5)
                    plt.text(0, res, "{0}:{1} {2:.2f}".format(i+1, i, res), color=c, horizontalalignment='right', zorder=5)
                for i in range(1, 4):
                    if i == 2: return
                    res = mmr(meana, (i+2)/i)
                    plt.axhline(y=res, linewidth=1, color=c, zorder=5)
                    plt.text(0, res, "{0}:{1} {2:.2f}".format(i+2, i, res), color=c, horizontalalignment='right', zorder=5)
                for i in range(1, 3):
                    if i == 3: return
                    res = mmr(meana, (i+3)/i)
                    plt.axhline(y=res, linewidth=1, color=c, zorder=5)
                    plt.text(0, res, "{0}:{1} {2:.2f}".format(i+3, i, res), color=c, horizontalalignment='right', zorder=5)

            plt.plot([], [], c=c, label="Planet {0}".format(ind))
        else:
            plt.plot([], [], c=c, label="Particle {0}".format(ind))

    do_for(plot_qQ)
    plt.title(sys.argv[1])
    plt.xlabel(timelabel)
    plt.ylabel("a, q, Q (AU)")
    plt.legend()

dt = times[1] - times[0]

if args["--plot-ftrect"]:
    fig, axes = plt.subplots(2, 1, sharex=True)

    def plot_hp(data, c, ind, is_planet):
        notnan = np.logical_not(np.isnan(data[0]))
        if not np.any(notnan):
            return

        taxis = np.fft.fftfreq(len(times[notnan]), dt / 365) * 360 * 3600
        half = len(times[notnan]) // 2

        axes[0].plot(taxis[1:half], np.abs(np.fft.fft(data[2, notnan] * np.sin(data[3, notnan]))[1:half]), c=c)
        axes[1].plot(taxis[1:half], np.abs(np.fft.fft(data[1, notnan] * np.sin(data[4, notnan] + data[3, notnan]))[1:half]), c=c)
        if is_planet:
            axes[0].plot([], [], c=c, label="Planet {0}".format(ind))
        else:
            axes[0].plot([], [], c=c, label="Particle {0}".format(ind))


    do_for(plot_hp)
    axes[0].legend()
    axes[0].set_title("p")
    axes[0].set_yscale("log")
    axes[0].set_xscale("log")
    axes[1].set_title("h")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("arcsec / yr")
    axes[0].legend()

if args["--plot-aei"]:
    fig, axes = plt.subplots(3, 1, sharex=True)

    def plot_aei(data, c, ind, is_planet):
        axes[0].plot(stimes, data[0] - data[0].mean(), c=c)
        axes[1].plot(stimes, data[1], c=c)
        axes[2].plot(stimes, data[2], c=c)
        if is_planet:
            axes[0].plot([], [], c=c, label="Planet {0}".format(ind))
        else:
            axes[0].plot([], [], c=c, label="Particle {0}".format(ind))

    axes[0].set_title("δa")
    axes[1].set_title("e")
    axes[2].set_title("i")
    axes[2].set_xlabel(timelabel)
    axes[0].legend()

if args["--plot-e-smooth"]:
    def moving_average(a, n=3):
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n
    fig, axes = plt.subplots(1, 1, sharex=True)

    def plot_aei(data, c, ind, is_planet):
        axes.plot(stimes[9:], moving_average(data[1], 10), c=c)

        if is_planet:
            axes.plot([], [], c=c, label="Particle {0}".format(ind))
        else:
            axes.plot([], [], c=c, label="Particle {0}".format(ind))

    axes.set_title("e")
    axes.set_xlabel(timelabel)
    axes.legend()

if args["--plot-angles"]:
    fig, axes = plt.subplots(3, 1, sharex=True)
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

    def plot_angles(data, c, ind, is_planet):
        axes[0].plot(stimes, normalize(data[3]), c=c)
        axes[1].plot(stimes, normalize(data[3]+data[4]), c=c)
        axes[2].plot(stimes, normalize(data[5]), c=c)
        if is_planet:
            axes[0].plot([], [], c=c, label="Planet {0}".format(ind))
        else:
            axes[0].plot([], [], c=c, label="Particle {0}".format(ind))

    axes[0].set_title("dΩ/dt (rad/yr)")
    axes[1].set_title("d(Ω+ω)/dt")
    axes[2].set_title("df/dt")
    axes[2].set_xlabel(timelabel)
    axes[0].legend()

plt.show()
