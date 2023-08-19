#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate Fig S15 a
"""
import numpy as np
from matplotlib import pylab as plt

path = 'data/'
c1 = path + "prion-data-model-tloss-lam1p75-uni-v1.txt"
c2 = path + "prion-data-model-tloss-lam1p75-uni-v1-small-seed.txt"
c3 = path + "prion-data-model-tloss-lam1p75-uni-v1-large-seed.txt"

data1 = np.loadtxt(c1)
data2 = np.loadtxt(c2)
data3 = np.loadtxt(c3)


timestep = 10
N = 200



count1, bins_count1 = np.histogram(data1, bins=200)
pdf1 = count1 / sum(count1)
cdf1 = np.cumsum(pdf1)


count2, bins_count2 = np.histogram(data2, bins=200)
pdf2 = count2 / sum(count2)
cdf2 = np.cumsum(pdf2)


count3, bins_count3 = np.histogram(data3, bins=200)
pdf3 = count3 / sum(count3)
cdf3 = np.cumsum(pdf3)


cmap=plt.get_cmap('Blues')  

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(bins_count3[1:], 1 - cdf3, color = 'C0', linewidth = 3, label=r'Initial seed size = $2s_{0}$')
ax.plot(bins_count1[1:], 1 - cdf1, '-',color = 'k', linewidth = 3,  label=r'Initial seed size = $s_{0}$ ')
ax.plot(bins_count2[1:], 1 - cdf2, '-', color = 'C3', linewidth = 3, label=r'Initial seed size = $s_{0}/2$')




ax.set_xlabel('Time (min)', fontsize = 15)
ax.set_ylabel('Fraction of cells with prions', fontsize = 15)
plt.title('Different initial seed size')
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
ax.set_xlim([0,500])
ax.set_ylim([0.01, 1])
ax.legend(loc = (0.5,0.55), framealpha=1)
plt.savefig('fig-SI-initial-seed.pdf')







