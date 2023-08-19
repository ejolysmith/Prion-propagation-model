#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates FigS19d
"""

import numpy as np
from matplotlib import pylab as plt
import statistics

path = 'data_time-rel-loss/'
c1 = path + "prion-data-model-6-time-rel-slow-stable-part1.txt"
#c1 = path + "prion-data-model-6-time-rel-slow-stable-v2.txt"
data = np.loadtxt(c1)
for i in range(2,51):
    print(i)
    c1 = path + "prion-data-model-6-time-rel-slow-stable-part" + str(i) + ".txt"
    data = np.concatenate((data, np.loadtxt(c1)))


ind = np.zeros(len(data))

for i in range(len(data)-1):
    ind[i] = data[:,2][i+1] - data[:,2][i]

indices_of_loss = np.where(ind == 1)[0] 
#indices_of_loss = indices_of_loss[0] - np.ones(len(indices_of_loss[0]))

    
indices_of_switch = (np.where(ind == -1))[0]

for j in range(int(indices_of_switch[0])):
    a = int(indices_of_switch[0])
    b = int(indices_of_loss[0])
    data[j,1] = data[j,1] - data[b,1]
    
for i in range(len(indices_of_switch)-1):
    a = int(indices_of_switch[i])
    b = int(indices_of_loss[i+1])
    B = data[b,1]
    for j in range(int(indices_of_switch[i+1] - indices_of_switch[i])):
        data[a,1] = data[a,1] - B
        a = a + 1
        #b = b - 1
        

number_of_bins = 9

time_bins1 = np.linspace(-400, 10, number_of_bins)
med_bins1 = []

for i in range(len(time_bins1)-1):
    li = []
    for j in range(len(data)):
        if time_bins1[i] <= data[j,1] < time_bins1[i+1] : 
            li.append(data[j,5])
    med_bins1.append(statistics.mean(li))
    print(len(li))
            

med_bins2 = []

for i in range(len(time_bins1)-1):
    li = []
    for j in range(len(data)):
        if time_bins1[i] <= data[j,1] < time_bins1[i+1] : 
            li.append(data[j,4])
    med_bins2.append(statistics.mean(li))
    print(len(li))
            
        


time_bins3 = np.linspace(-400, 10, number_of_bins)
med_bins3 = []

for i in range(len(time_bins3)-1):
    med_bins3.append(1)
    print(len(li))

mean1 = 1#med_bins1[0]
mean2 = 1
mean3 = 1
plt.plot(np.array(time_bins1[1:] - time_bins1[1]/2 + time_bins1[0]/2), np.array(med_bins1)/np.array(med_bins2), linewidth = 5, color = 'b',label='Average chaperone abundance = 1')
#plt.plot(np.array(time_bins2[1:] - time_bins2[1]/2 + time_bins2[0]/2), np.array(med_bins2)/mean2, linewidth = 5, color = 'r',label='Average chaperone abundance = 10')
#plt.plot(np.array(time_bins1[1:] - time_bins1[1]/2 + time_bins1[0]/2), np.array(med_bins3)/mean3, linewidth = 5, color = 'k',label='Average chaperone abundance = 100')

plt.xlabel("Time relative to loss (min)", fontsize = 16)
#plt.yticks([0,0.5,1,1.5,2], fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

plt.ylabel("Mean aggregate size", fontsize = 16)
#plt.ylim([0.5,1.5])
plt.legend()

plt.savefig('size-time-rel-slow-stable-v1.pdf')

#plt.title("alpha = 1e-3 and gamma = 3e-3")
#plt.xlim([-80, 10])
