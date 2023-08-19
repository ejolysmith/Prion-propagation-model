#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates Fig S18a

"""

import numpy as np
from matplotlib import pylab as plt


path = 'data_single-chaperones-traces/'
c1 = path + "prion-data-single-chaperones-slow-long-v1.txt"


data1 = np.loadtxt(c1)
fig, ax = plt.subplots(1, 1)
t = data1[:,0]
c = data1[:,1]
ax.plot(t,c, label = r'Slow pulses')
ax.set_ylabel("Chaperone  conc.  \n (arb. units)", fontsize = 16)
ax.set_xlabel("Time (min)", fontsize = 16)
plt.xticks(fontsize = 14)
plt.legend()
plt.yticks([0,25,59,75,100,125], fontsize = 14)
#plt.yticks([0,1,2,3], ["0.0", "1.0", "2.0", "3.0"], fontsize = 14)
plt.tight_layout()
plt.savefig('single-1.pdf')



fig, ax = plt.subplots(1, 1)
t = data1[:,0]
c = np.ones(len(t))*70
ax.plot(t,c, label = r'Slow pulses')
ax.set_ylabel("Chaperone  conc.  \n (arb. units)", fontsize = 16)
ax.set_xlabel("Time (min)", fontsize = 16)
plt.xticks(fontsize = 14)
plt.legend()
plt.yticks([0,25,59,75,100,125], fontsize = 14)
plt.tight_layout()
plt.savefig('single-2.pdf')




