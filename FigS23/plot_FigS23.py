"""
generates FigS23

"""

import numpy as np
from matplotlib import pylab as plt

c1 = "prion-simulation-data-model-6-divisions-relative-to-loss.txt"

data = np.loadtxt(c1)


    
    
fig = plt.figure()
plt.plot(data[:,-3], data[:,-1])
plt.xlabel("Time relative to loss (min)", fontsize=16)
plt.ylabel("Deterministic" + "\n" + " aggregate abundance", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig("fig-deterministic-abundance.pdf")


fig = plt.figure()
plt.plot(data[:,-3], data[:,-2])
plt.xlabel("Time relative to loss (min)", fontsize=16)
plt.ylabel("Deterministic" + "\n" + " aggregate concentration (a.u.)", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig("fig-deterministic-conc.pdf")

