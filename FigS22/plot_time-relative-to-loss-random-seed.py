"""
Generates top left plot in FigS22

"""

import numpy as np
from matplotlib import pylab as plt
import statistics

path = 'data_time-rel-to-loss/'
c1 = path + "prion-data-model-6-time-rel-v1.txt"

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
        




number_of_bins = 20

time_bins1 = np.linspace(-80, 10, number_of_bins)
med_bins1 = []

for i in range(len(time_bins1)-1):
    li = []
    for j in range(len(data)):
        if time_bins1[i] <= data[j,1] < time_bins1[i+1] : 
            li.append(data[j,0])
    med_bins1.append(statistics.median(li))
    print(len(li))
            


c1 = path +"prion-data-model-6-time-rel-init-v1.txt"

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
        




number_of_bins = 20

time_bins2 = np.linspace(-80, 10, number_of_bins)
med_bins2 = []

for i in range(len(time_bins2)-1):
    li = []
    for j in range(len(data)):
        if time_bins2[i] <= data[j,1] < time_bins2[i+1] : 
            li.append(data[j,0])
    med_bins2.append(statistics.median(li))
    print(len(li))
    
plt.plot(np.array(time_bins1[1:] - time_bins1[1]/2 + time_bins1[0]/2), med_bins1, linewidth = 5, color = 'k')
plt.plot(np.array(time_bins2[1:] - time_bins2[1]/2 + time_bins2[0]/2), med_bins2, linewidth = 5, color = 'darkorange')
plt.xlabel("Time relative to loss (min)", fontsize = 16)
plt.yticks([0,50,100,150,200], fontsize = 14)
plt.xticks(fontsize = 14)
plt.ylabel("Median protein" + "\n" +'concentration (a.u.)', fontsize = 16)
plt.ylim([0,200])
#plt.title("alpha = 1e-3 and gamma = 3e-3")
plt.savefig("fig-time-rel-to-loss-random-seed.pdf")
