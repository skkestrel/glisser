"""
============
3D animation
============

A simple example of an animated plot... In 3D!
"""
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import sys
import matplotlib.animation as animation

start = 0
interval = 20

if len(sys.argv) > 2:
    start = int(sys.argv[2])
if len(sys.argv) > 3:
    interval = int(sys.argv[3])

with open(sys.argv[1]) as f:
    x = []
    xp = []
    for line in f:
        vals = line.strip().split()

        xp.append((float(vals[0]), float(vals[1]), float(vals[2])))

data = [np.array(x).T, np.array(xp).T]

data = [data[1][:, ::10]]

def update_lines(num, dataLines, lines, txt):
    for line, data in zip(lines, dataLines):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, max(num-10, 0):num])
        line.set_3d_properties(data[2, max(num-10, 0):num])
    txt.set_text("t = {0}".format(num))
    return lines

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

# NOTE: Can't pass empty arrays into 3d version of plot()
lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1], label="")[0] for dat in data]

# Setting the axes properties
ax.set_xlim3d([-10.0, 10.0])
ax.set_xlabel('X')

ax.set_ylim3d([-10.0, 10.0])
ax.set_ylabel('Y')

ax.set_zlim3d([-0.5, 0.5])
ax.set_zlabel('Z')

ax.set_title('3D Test')

txt = ax.text(0, 0, 0.3, "")

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, update_lines, range(start, data[0].shape[1]), fargs=(data, lines, txt),
                                   interval=interval, blit=False)

plt.show()
