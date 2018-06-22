import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import struct
import matplotlib.style as style

style.use('ggplot')

times = []
planets = [[] for i in range(8)]
particlewatch = [int(sys.argv[i + 2]) for i in range(len(sys.argv) - 2)]
particles = [[] for i in range(len(particlewatch) * 2)]

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

with open(sys.argv[1], 'rb') as f:
    read = f.read(16)
    while len(read) == 16:
        time, npl = struct.unpack('=dQ', read)
        times.append(time)
        for i in range(npl):
            a, e = struct.unpack('=I3f', f.read(16))[1:3]
            planets[2*i].append(a)
            planets[2*i+1].append(e)
        npa, = struct.unpack('=Q', f.read(8))

        for index, partnum in enumerate(particlewatch):
            part = bsearch(f, npa, partnum)
            if part:
                particles[2*index].append(part[1])
                particles[2*index+1].append(part[2])
            else:
                particles[2*index].append(math.nan)
                particles[2*index+1].append(math.nan)

        f.seek(npa * 16, 1)
        read = f.read(16)

if times[-1] > 365e6:
    times = np.array(times) / 365e6
    plt.xlabel("Time (Myr)")
elif times[-1] > 365e3:
    times = np.array(times) / 365e3
    plt.xlabel("Time (kyr)")
else:
    times = np.array(times)
    plt.xlabel("Time (yr)")

planets = np.array(planets)
particles = np.array(particles)
cc = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i in range(4):
    c = cc[(i) % len(cc)]
    plt.plot(times, planets[2*i], c=c)
    plt.plot(times, planets[2*i] * (1. - planets[2*i+1] ), c=c)
    plt.plot(times, planets[2*i] * (1. + planets[2*i+1] ), c=c)
    plt.plot([], [], c=cc[i%len(cc)], label="Planet {0}".format(i+1))
for i in range(len(particlewatch)):
    c = cc[(i+4) % len(cc)]
    nnull = particles[2*i] != math.nan

    plt.plot(times[nnull], particles[2*i][nnull], c=c)
    plt.plot(times[nnull], (particles[2*i] * ( 1. - particles[2*i+1]) )[nnull], c=c)
    plt.plot(times[nnull], (particles[2*i] * ( 1. + particles[2*i+1]) )[nnull], c=c)
    plt.plot([], [], c=c, label="Particle {0}".format(particlewatch[i]))

plt.title(sys.argv[1])
plt.legend()
plt.show()
