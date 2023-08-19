#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: euan

This code takes in the simulated cell cycle position at time of loss data from 
the compiled prion-model_cell-cycle-position.c 
"""

import numpy as np
from matplotlib import pylab as plt

c1 = "prion-data-model-cell-cycle-time-lam1p75-uni-v1.txt"

data = np.loadtxt(c1)

bins = np.linspace(0,1,11)

binned_data = np.zeros(10)

for i in range(10):
    for j in range(len(data)):
        if bins[i] <= data[j] < bins[i+1]:
            binned_data[i] = binned_data[i] + 1
            
            
plt.hist(data, bins = bins , weights=np.ones(len(data)) / len(data))
plt.xlabel('Cell cycle position at time of loss', fontsize = 18)
plt.ylabel('Fraction of simulations', fontsize = 18)
plt.xticks(fontsize = 17)
plt.yticks(fontsize = 17)

            
        