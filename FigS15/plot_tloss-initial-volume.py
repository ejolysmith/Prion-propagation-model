#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates Fig S15b

"""
import numpy as np
from matplotlib import pylab as plt

path = 'data/'
c1 = path + "prion-data-model-tloss-lam1p75-uni-v1.txt"
c2 = path + "prion-data-model-tloss-lam1p75-exp-v1.txt"
c3 = path + "prion-data-model-tloss-lam1p75-iden-v1.txt"

data1 = np.loadtxt(c1)
data2 = np.loadtxt(c2)
data3 = np.loadtxt(c3)


timestep = 2
N = 200

count1, bins_count1 = np.histogram(data1, bins=200)
pdf1 = count1 / sum(count1)
cdf1 = np.cumsum(pdf1)


count2, bins_count2 = np.histogram(data2, bins=200)
pdf2 = count2 / sum(count2)
cdf2 = np.cumsum(pdf2)


count3, bins_count3 = np.histogram(data3, bins=1000)
pdf3 = count3 / sum(count3)
cdf3 = np.cumsum(pdf3)


cmap=plt.get_cmap('Blues')  

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
data3 = np.append(data3,0)  #Added to close the curve in the plot. This is because in this case (synchronized divisions)
#there are no losses prior to the first division.

ax.plot(bins_count3[1:], 1 - cdf3, color = 'C0', linewidth = 3, label=r'Synchronized cells $V(0) = V_{0}$') 
ax.plot(bins_count1[1:], 1 - cdf1, '-',color = 'k', linewidth = 3, label=r'Uniform distribution ')
ax.plot(bins_count2[1:], 1 - cdf2, '-', color = 'C3', linewidth = 3, label=r'Exponential distribution$')



ax.set_xlabel('Time (min)', fontsize = 15)
ax.set_ylabel('Fraction of cells with prions', fontsize = 15)
plt.title('Different initial cell size distributions')
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
ax.set_xlim([0,500])
ax.set_ylim([0.01, 1])
ax.legend(loc = (0.35,0.55), framealpha=1)

plt.savefig('fig-SI-initial-size.pdf')
