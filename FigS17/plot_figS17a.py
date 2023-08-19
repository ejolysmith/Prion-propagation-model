#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates fig S17 a

"""

import numpy as np
from matplotlib import pylab as plt


path = 'data_single-chaperones-traces/'
c1 = path + "prion-data-single-chaperones-1-v1.txt"
c2 = path + "prion-data-single-chaperones-10-v1.txt"
c3 = path + "prion-data-single-chaperones-100-v1.txt"


data1 = np.loadtxt(c1)

fig, ax = plt.subplots(1, 1, figsize = (5,3))
t = data1[:,0]
c = data1[:,1]
ax.plot(t,c, label = r'Average chaperone abundance = 1')
ax.set_ylabel("Chaperone  conc.  \n (arb. units)", fontsize = 16)
ax.set_xlabel("Time (min)", fontsize = 16)
plt.xticks(fontsize = 14)
plt.legend()
plt.yticks([0,1,2,3], ["0.0", "1.0", "2.0", "3.0"], fontsize = 14)
plt.tight_layout()
#plt.savefig('single-1.pdf')


data2 = np.loadtxt(c2)

fig, ax = plt.subplots(1, 1, figsize = (5,3))
t = data2[:,0]
c = data2[:,1]
ax.plot(t,c, label = r'Average chaperone abundance = 10')
ax.set_ylabel("Chaperone  conc.  \n (arb. units)", fontsize = 16)
ax.set_xlabel("Time (min)", fontsize = 16)
ax.set_ylim([0,15])
plt.legend()
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.yticks([0,5,10,15], ["0.0", "5.0", "10", "15"], fontsize = 14)

plt.tight_layout()
#plt.savefig('single-2.pdf')


data3 = np.loadtxt(c3)

fig, ax = plt.subplots(1, 1, figsize = (5,3))
t = data3[:,0]
c = data3[:,1]
ax.plot(t,c, label = r'Average chaperone abundance = 100')
ax.set_ylabel("Chaperone  conc. \n (arb. units)", fontsize = 16)
ax.set_xlabel("Time (min)", fontsize = 16)
ax.set_ylim([0,90])
plt.legend()
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.yticks([0,20,40,60,80], ["0.0", "20", "40", "60", "80"], fontsize = 14)

plt.tight_layout()
#plt.savefig('single-3.pdf')
