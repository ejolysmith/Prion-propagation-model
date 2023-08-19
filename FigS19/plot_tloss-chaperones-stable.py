#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates Fig S19bc

"""

import numpy as np
from matplotlib import pylab as plt


path = 'data_tloss/'
c1 = path + "prion-data-model-6-tloss-v1-slow-very-stable.txt"

data1 = np.loadtxt(c1)


timestep = 10
N = 100



count1, bins_count1 = np.histogram(data1, bins=500)
pdf1 = count1 / sum(count1)
cdf1 = np.cumsum(pdf1)



cmap=plt.get_cmap('Blues')  

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(bins_count1[1:], 1 - cdf1, color = cmap(0.5))



ax.set_xlabel('Time (min)', fontsize = 16)
ax.set_ylabel('Fraction of cells with prions', fontsize = 16)
ax.set_xlim([0,10000])
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

#ax.set_ylim([1e-3, 1])
plt.savefig( "fig-tloss-chaperones-stable.pdf" )




    


