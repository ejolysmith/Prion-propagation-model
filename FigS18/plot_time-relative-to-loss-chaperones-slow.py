"""
Generates Fig S18b

"""

import numpy as np
from matplotlib import pylab as plt
import statistics


path = 'data_time-rel-loss/'
c1 = path + "prion-data-model-6-time-rel-slow-v1-part1.txt"
data = np.loadtxt(c1)
for i in range(2,31):
    print(i)
    c1 = path + "prion-data-model-6-time-rel-slow-v1-part" + str(i) + ".txt"
    data = np.concatenate((data, np.loadtxt(c1)))

ind = np.zeros(len(data))

for i in range(len(data)-1):
    ind[i] = data[:,2][i+1] - data[:,2][i]

indices_of_loss = np.where(ind == 1)[0] 

    
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
        

number_of_bins = 8

time_bins1 = np.linspace(-400, 10, number_of_bins)
med_bins1 = []

for i in range(len(time_bins1)-1):
    li = []
    for j in range(len(data)):
        if time_bins1[i] <= data[j,1] < time_bins1[i+1] : 
            li.append(data[j,3])
    med_bins1.append(statistics.mean(li))
    print(len(li))
            
    


mean1 = med_bins1[0]
plt.plot(np.array(time_bins1[1:] - time_bins1[1]/2 + time_bins1[0]/2), np.array(med_bins1)/mean1, linewidth = 5, color = 'b',label='Average chaperone abundance = 1')
plt.xlabel("Time relative to loss (min)", fontsize = 16)
plt.yticks([0,0.5,1,1.5,2], fontsize = 14)
plt.xticks(fontsize = 14)
plt.ylabel("Normalized mean chaperone" + "\n" +'concentration (a.u.)', fontsize = 16)
plt.ylim([0,2])
#plt.xticks([-400, -300, -200, -100, 0])
plt.legend()
plt.savefig("fig-time-relative-to-loss-chaperones-slow.pdf")
