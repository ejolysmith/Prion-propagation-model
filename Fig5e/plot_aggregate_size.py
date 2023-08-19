#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: euan

This code takes in the simulated data from 
the compiled prion_model_dist.c .  This data is the distribution shown in 
Fig. S9b (discussed in caption), but it needs to be normalized (which is done here). 

"""


import numpy as np
from matplotlib import pylab as plt


path = 'data/'


#Data obtained from prion-model_dist.c by varying gamma and alpha 
#with the following respective pairs (gamma, alpha):
#(1.12200000e-03, 1.97000000e-04)
#(0.001902, 0.0009097)
#(0.001805, 0.0008194)
#(0.001707, 0.000729)
#(0.001512, 0.0005484)
c1 = path + "prion-data-model-dist-lam1p75-uni-v1.txt"
c2 = path + "prion-data-model-dist-lam1p75-uni-v2.txt"
c3 = path + "prion-data-model-dist-lam1p75-uni-v3.txt"
c4 = path + "prion-data-model-dist-lam1p75-uni-v4.txt"
c5 = path + "prion-data-model-dist-lam1p75-uni-v5.txt"

data1 = np.loadtxt(c1)
data2 = np.loadtxt(c2)
data3 = np.loadtxt(c3)
data4 = np.loadtxt(c4)
data5 = np.loadtxt(c5)

timestep = 10
N = 100 


for i in range(len(data1[:,0])):
    data1[:,2][i] = data1[:,2][i]/data1[:,1][i]/10
    sum = np.sum(data1[i,:])
    for j in range(3, 2002):
        if data1[i,1] > 0 :
            data1[i,j] = data1[i,j]/sum
        

for i in range(len(data2[:,0])):
    data2[:,2][i] = data2[:,2][i]/data2[:,1][i]/10
    sum = np.sum(data2[i,:])
    for j in range(3, 2002):
        if data2[i,1] > 0 :
            data2[i,j] = data2[i,j]/sum
            
for i in range(len(data3[:,0])):
    data3[:,2][i] = data3[:,2][i]/data3[:,1][i]/10
    sum = np.sum(data3[i,:])
    for j in range(3, 2002):
        if data3[i,1] > 0 :
            data3[i,j] = data3[i,j]/sum
            
for i in range(len(data4[:,0])):
    data4[:,2][i] = data4[:,2][i]/data4[:,1][i]/10
    sum = np.sum(data4[i,:])
    for j in range(3, 2002):
        if data4[i,1] > 0 :
            data4[i,j] = data4[i,j]/sum
            
for i in range(len(data5[:,0])):
    data5[:,2][i] = data5[:,2][i]/data5[:,1][i]/10
    sum = np.sum(data5[i,:])
    for j in range(3, 2002):
        if data5[i,1] > 0 :
            data5[i,j] = data5[i,j]/sum




fig = plt.figure()



plt.plot(range(1, 2001), data5[-1,3:], linewidth = 5)
plt.plot(range(1, 2001), data4[-1,3:], linewidth = 5)
plt.plot(range(1, 2001), data3[-1,3:], linewidth = 5)
plt.plot(range(1, 2001), data2[-1,3:], linewidth = 5)
plt.plot(range(1, 2001), data1[-1,3:], linewidth = 5)


plt.xlabel("Aggregate size", fontsize = 15)
plt.ylabel("Fraction of total aggregates", fontsize = 15)
plt.xlim([0, 80])
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
#plt.savefig("fig-dist-rel.pdf")