"""
@author: euan

This code takes in the data from 
the compiled prion-time_relative_to_loss.c. This data is a long array of
protein concentration measurements taken at every time a reaction occurs over 
1000 simulations. This code orders these measurements according to their time 
relative to the loss event (in bins), and plots the average of concentration 
relative to loss. 

"""

import numpy as np
from matplotlib import pylab as plt
import statistics

path = 'data/'
c1 = path + "prion-data-model-time-rel-lam1p75-uni-v1.txt"
data = np.loadtxt(c1)

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
        


number_of_bins = 30

time_bins = np.linspace(-80, 10, number_of_bins)
med_bins = []

for i in range(len(time_bins)-1):
    li = []
    for j in range(len(data)):
        if time_bins[i] <= data[j,1] < time_bins[i+1] : 
            li.append(data[j,0])
    med_bins.append(statistics.median(li))
    print(len(li))
            

plt.plot(np.array(time_bins[1:] - time_bins[1]/2 + time_bins[0]/2), med_bins, linewidth = 5, color = 'darkorange')
plt.xlabel("Time relative to loss (min)", fontsize = 25)
plt.yticks([0,50,100,150,200], fontsize = 20)
plt.xticks(fontsize = 20)
plt.ylabel("Median protein" + "\n" +'concentration (a.u.)', fontsize = 25)
plt.ylim([0,200])
#plt.title("alpha = 1e-3 and gamma = 3e-3")
#plt.savefig("fig-time-rel.pdf")
