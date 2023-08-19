#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates Fig S16
@author: euan
"""
from __future__ import division



import numpy as np
from matplotlib import pylab as plt


c1 = "data_seed2.txt"


data = np.loadtxt(c1)

        


# Generate fake data
x = data[:,7]
y = data[:,6]

# Calculate the point density
xy = np.vstack([x,y])


z = data[:,0]

fig, ax = plt.subplots()
density = ax.scatter(x, y, c=z, s=100, vmin = 0.0, vmax = 500, cmap="Blues")#0.125)
#density2 = ax.scatter([1.3*4e-4], [0.8*2e-4], c=[0], s=100, vmin = 0.0, vmax = 1000)#0.125)

plt.xscale('log')
plt.yscale('log')
#plt.ylim([1e-4, 1e-2])
#plt.xlim([1e-4, 1e-2])
cb = fig.colorbar(density, label='Average time of loss')    
plt.ylim([2*1e-5, 1e-2])
plt.xlim([2*1e-5, 1e-1])
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.xlabel(r'Elongation rate $\gamma V_{0}$', fontsize = 15)

cb.set_label(label='Average time of loss (min)',size = 15, rotation=-90, labelpad = 17)
cb.set_ticks([0,100,200,300,400,500])
cb.set_ticklabels([0, 100, 200, 300, 400, '>500'])

plt.xlabel(r"Elongation rate $\gamma V_{0}$ (min$^{-1}$)", fontsize = 15)
plt.ylabel(r"Fragmentation rate $\alpha V_{0}$ (min$^{-1}$)", fontsize = 15)
alphadot = 1*2.48000000e-04
gammadot = 1.08300000e-03
plt.scatter([gammadot], [alphadot], color = 'darkorange', s = 200)


plt.savefig('fig-SI-density-time-of-loss-seed2.pdf')

xs = np.linspace(1e-4, 1e-1, 1000)
#plt.plot(xs, func(xs))

z = data[:,3]

fig, ax = plt.subplots()
density = ax.scatter(x, y, c=z, s=100, vmin = 0.0, vmax = 0.5, cmap="Oranges")#0.125)
#density2 = ax.scatter([1.3*4e-4], [0.8*2e-4], c=[0], s=100, vmin = 0.0, vmax = 0.5)#0.125)

plt.xscale('log')
plt.yscale('log')
plt.ylim([2*1e-5, 1e-2])
plt.xlim([2*1e-5, 1e-1])
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
cb = fig.colorbar(density, label='Absolute partitioning error')    
cb.set_label(label='Absolute partitioning error',size = 15, rotation=-90, labelpad = 17)
plt.xlabel(r'Elongation rate $\gamma V_{0}$', fontsize = 15)

cb.set_label(label='Absolute partitioning error',size = 15, rotation=-90, labelpad = 17)
cb.set_ticks([0,0.1,0.2,0.3,0.4,0.5])
cb.set_ticklabels([0, 0.1, 0.2, 0.3, 0.4, '>0.5'])

plt.xlabel(r"Elongation rate $\gamma V_{0}$ (min$^{-1}$)", fontsize = 15)
plt.ylabel(r"Fragmentation rate $\alpha V_{0}$ (min$^{-1}$)", fontsize = 15)
alphadot = 1*2.48000000e-04
gammadot = 1.08300000e-03
plt.scatter([gammadot], [alphadot], color = 'darkorange', s = 200)

plt.savefig('fig-SI-density-abs-part-error-seed2.pdf')



z = data[:,1]

fig, ax = plt.subplots()
density = ax.scatter(x, y, c=z, s=100, vmin = 0.0, vmax = 10, cmap="Purples")#0.125)
#density2 = ax.scatter([1.3*4e-4], [0.8*2e-4], c=[0], s=100, vmin = 0.0, vmax = 0.5)#0.125)

plt.xscale('log')
plt.yscale('log')
plt.ylim([2*1e-5, 1e-2])
plt.xlim([2*1e-5, 1e-1])
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
cb = fig.colorbar(density, label='Average number of aggregates')    
cb.set_label(label='Average number of aggregates',size = 15, rotation=-90, labelpad = 17)
plt.xlabel(r'Elongation rate $\gamma V_{0}$', fontsize = 15)

cb.set_label(label='Average number of aggregates',size = 15, rotation=-90, labelpad = 17)
cb.set_ticks([0, 2, 4, 6, 8, 10])
cb.set_ticklabels([0, 2, 4, 6, 8, '>10'])

plt.xlabel(r"Elongation rate $\gamma V_{0}$ (min$^{-1}$)", fontsize = 15)
plt.ylabel(r"Fragmentation rate $\alpha V_{0}$ (min$^{-1}$)", fontsize = 15)
alphadot = 1*2.48000000e-04
gammadot = 1.08300000e-03
plt.scatter([gammadot], [alphadot], color = 'darkorange', s = 200)

plt.savefig('fig-SI-density-number-of_aggregates-seed2.pdf')


z = data[:,2]

fig, ax = plt.subplots()
density = ax.scatter(x, y, c=z, s=100, vmin = 0.0, vmax = 80, cmap="Greens",)#0.125)
#density2 = ax.scatter([1.3*4e-4], [0.8*2e-4], c=[0], s=100, vmin = 0.0, vmax = 0.5)#0.125)

plt.xscale('log')
plt.yscale('log')
plt.ylim([2*1e-5, 1e-2])
plt.xlim([2*1e-5, 1e-1])
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
cb = fig.colorbar(density, label='Average aggregate size')    
cb.set_label(label='Average aggregate size',size = 15, rotation=-90, labelpad = 17)
plt.xlabel(r'Elongation rate $\gamma V_{0}$', fontsize = 15)

cb.set_label(label='Average aggregate size',size = 15, rotation=-90, labelpad = 17)
cb.set_ticks([0, 20, 40, 60, 80])
cb.set_ticklabels([0, 20, 40, 60, '>80'])

plt.xlabel(r"Elongation rate $\gamma V_{0}$ (min$^{-1}$)", fontsize = 15)
plt.ylabel(r"Fragmentation rate $\alpha V_{0}$ (min$^{-1}$)", fontsize = 15)
alphadot = 1*2.48000000e-04
gammadot = 1.08300000e-03
plt.scatter([gammadot], [alphadot], color = 'darkorange', s = 200)

plt.savefig('fig-SI-density-aggregate-size-seed2.pdf')


