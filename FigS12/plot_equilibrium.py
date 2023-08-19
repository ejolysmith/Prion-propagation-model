#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
To plot Fig.S12
"""

import numpy as np
from matplotlib import pylab as plt

c1 = "prion-data-model-equilibrium-v1.txt"

data = np.loadtxt(c1)
plt.plot(data[:,0]/60, data[:,1], linewidth = 4)
#370 is taken from the printed output of prion-model_equilibrium.c
plt.plot(np.ones(100)*370./60, np.linspace(0,1000,100),  '--', color = 'k', linewidth = 1)
plt.ylim([0, 1000])
plt.xlim([0,15])
plt.yticks(fontsize = 12)
plt.xticks([0,5,10,15], fontsize = 12)
plt.xlabel('Time (h)', fontsize = 15)
plt.ylabel('Mean protein \n concentration (arb. units)', fontsize = 15)
plt.tight_layout()
plt.savefig('fig-SI-equilibrium.pdf')

