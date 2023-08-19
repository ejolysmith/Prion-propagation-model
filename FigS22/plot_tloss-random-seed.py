#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates top right plot in FigS22
"""
import numpy as np
from matplotlib import pylab as plt



path = 'data_tloss/'
c1 = path + "prion-data-model-6-tloss-v1.txt"
c2 = path + "prion-data-model-6-tloss-init-v1.txt"

data1 = np.loadtxt(c1)
data2 = np.loadtxt(c2)

timestep = 10
N = 100



count1, bins_count1 = np.histogram(data1, bins=100)
pdf1 = count1 / sum(count1)
cdf1 = np.cumsum(pdf1)


count2, bins_count2 = np.histogram(data2, bins=100)
pdf2 = count2 / sum(count2)
cdf2 = np.cumsum(pdf2)



cmap=plt.get_cmap('Blues')  

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#ax.set_yscale('log')
#ax.plot(bins_count[1:], pdf, color="red", label="PDF")
ax.plot(bins_count1[1:], 1 - cdf1, linewidth = 2,color = 'darkorange', label='Constant initial seed size')
ax.plot(bins_count2[1:], 1 - cdf2, linewidth = 2,color = 'k', label='Random initial seed size')



ax.set_xlabel('Time (min)', fontsize = 16)
ax.set_ylabel('Fraction of cells with prions', fontsize = 16)
ax.set_xlim([0,700])
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

#ax.set_ylim([1e-3, 1])
plt.legend()
plt.savefig( "fig-tloss-random-seed.pdf" )


    


