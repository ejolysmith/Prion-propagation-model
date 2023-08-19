#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: euan

This code takes in the simulated time-of-loss data from 
the compiled prion-model_tloss.c 

"""
import numpy as np
from matplotlib import pylab as plt


#Data obtained from prion-model_tloss.c by varying gamma and alpha 
#with the following respective pairs (gamma, alpha):
#(1.12200000e-03, 1.97000000e-04)
#(0.001902, 0.0009097)
#(0.001805, 0.0008194)
#(0.001707, 0.000729)
#(0.001512, 0.0005484)


path = 'data/'
c1 = path + "prion-data-model-tloss-lam1p75-uni-v1.txt"
c2 = path + "prion-data-model-tloss-lam1p75-uni-v2.txt"
c3 = path + "prion-data-model-tloss-lam1p75-uni-v3.txt"
c4 = path + "prion-data-model-tloss-lam1p75-uni-v4.txt"
c5 = path + "prion-data-model-tloss-lam1p75-uni-v5.txt"


data1 = np.loadtxt(c1)
data2 = np.loadtxt(c2)
data3 = np.loadtxt(c3)
data4 = np.loadtxt(c4)
data5 = np.loadtxt(c5)


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


count4, bins_count4 = np.histogram(data4, bins=200)#
pdf4 = count4 / sum(count4)
cdf4 = np.cumsum(pdf4)


count5, bins_count5 = np.histogram(data5, bins=200)
pdf5 = count5 / sum(count5)
cdf5 = np.cumsum(pdf5)


cmap=plt.get_cmap('Blues')  

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(bins_count5[1:], 1 - cdf5, color = cmap(0.7), linewidth = 3, label=r' ')
ax.plot(bins_count4[1:], 1 - cdf4, color = cmap(0.6), linewidth = 3, label=r' ')
ax.plot(bins_count3[1:], 1 - cdf3, color = cmap(0.5), linewidth = 3, label=r' ')
ax.plot(bins_count2[1:], 1 - cdf2, color = cmap(0.4), linewidth = 3, label=r' ')
ax.plot(bins_count1[1:], 1 - cdf1, color = cmap(0.3), linewidth = 3, label=r' ')



ax.set_xlabel('Time (min)', fontsize = 15)
ax.set_ylabel('Fraction of cells with prions', fontsize = 15)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
ax.set_xlim([0,1000])
ax.set_ylim([0.00, 1])
ax.legend(loc = (0.5,0.55), framealpha=1)
#plt.savefig('figure-tloss.pdf')



    


