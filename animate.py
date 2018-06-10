"""
============
3D animation
============

A simple example of an animated plot... In 3D!
"""
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

with open("animation.out", "r") as f:
	arrs = [[] for _ in range(505)]
	for index, line in enumerate(f):
		arrs[index % 505].append([float(x) for x in line.strip().split()])

arrs = [np.array(x).T for x in arrs]
		
def update_lines(num, dataLines, lines):
    for i, (line, data) in enumerate(zip(lines, dataLines)):
        n = 20 if i < 5 else 2
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, max(num-n, 0):num])
        line.set_3d_properties(data[2, max(num-n, 0):num])
    return lines

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in arrs]
# Creating fifty line objects.
# NOTE: Can't pass empty arrays into 3d version of plot()

# Setting the axes properties
ax.set_xlim3d([-30.0, 30.0])
ax.set_xlabel('X')

ax.set_ylim3d([-30.0, 30.0])
ax.set_ylabel('Y')

ax.set_zlim3d([-30.0, 30.0])
ax.set_zlabel('Z')

ax.set_title('3D Test')

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, update_lines, 500, fargs=(arrs, lines),
                                   interval=50, blit=False)

plt.show()
