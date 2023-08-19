#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates Fig S17c
"""
import numpy as np
from matplotlib import pylab as plt


path = 'data_tloss/'

c1 = path + "prion-data-model-6-tloss-1-v1.txt"
c2 = path + "prion-data-model-6-tloss-10-v1.txt"
c3 = path + "prion-data-model-6-tloss-100-v1.txt"



data1 = np.loadtxt(c1)
data2 = np.loadtxt(c2)
data3 = np.loadtxt(c3)



timestep = 10
N = 100
count1, bins_count1 = np.histogram(data1, bins=100)
pdf1 = count1 / sum(count1)
cdf1 = np.cumsum(pdf1)


count2, bins_count2 = np.histogram(data2, bins=100)
pdf2 = count2 / sum(count2)
cdf2 = np.cumsum(pdf2)


count3, bins_count3 = np.histogram(data3, bins=100)
pdf3 = count3 / sum(count3)
cdf3 = np.cumsum(pdf3)


cmap=plt.get_cmap('Blues')  

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(bins_count1[1:], 1 - cdf1, color = cmap(0.5), label='Average chaperone abundance = 1')
ax.plot(bins_count2[1:], 1 - cdf2, color = 'r', label='Average chaperone abundance = 10')
ax.plot(bins_count3[1:], 1 - cdf3, color = cmap(0.8), label='Average chaperone abundance = 100')



ax.set_xlabel('Time (min)', fontsize = 16)
ax.set_ylabel('Fraction of cells in with prions', fontsize = 16)
ax.set_xlim([0,700])
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

plt.legend()
plt.savefig("fig-tloss-chaperones.pdf")



    


