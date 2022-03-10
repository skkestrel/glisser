import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("temp_log.txt")
t, energy1, energy2 = data[:,0], data[:,1], data[:,2]

fig, ax = plt.subplots(figsize=(9,6))

ax.scatter(t, energy1, label='-GM/2a',s=0.2)
ax.scatter(t, energy2, label='v^2/2-GM/r',s=0.2)
ax.set_xlabel("Time (Days)")    
ax.set_ylabel("energy")  
ax.legend()
plt.savefig("energyyyy_4timestep.png",dpi=300)