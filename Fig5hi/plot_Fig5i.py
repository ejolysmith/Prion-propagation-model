"""
@author: euan

This code takes in the data from 
the compiled prion-divisions_relative_to_loss.c. This data is a long array of
partitioning errorstaken at every time a division occurs over 
10000 simulations. This code orders these divisions according relative to the 
loss event, and plots the average of partititioning error relative to loss. 

"""

import numpy as np
from matplotlib import pylab as plt

c1 = "data/prion-simulation-data-model-divisions-relative-to-loss.txt"

data = np.loadtxt(c1)


indices_of_loss = np.where(data[:,1] == 1)
indices_of_loss = (indices_of_loss[0] - np.ones(len(indices_of_loss[0])))

indices_of_switch = (np.where(data[:,1] == 10))[0]


for i in range(len(indices_of_loss) - 1):
    a = int(indices_of_loss[i+1])
    b = 0
    for j in range(int(indices_of_loss[i+1] - indices_of_switch[i])):
        data[a,1] = b
        #print(j)
        a = a - 1
        b = b - 1
        
     
divs = [-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,7,8,9,10]

part_error = np.zeros(len(divs))
abs_part_error = np.zeros(len(divs))

for i in range(len(divs)):
    where_i = np.where(data[:,1] == divs[i])[0]
    number_of_points = 0    
    for j in range(len(where_i)):
        index = int(where_i[j])
        part_error[i] = part_error[i] + data[index,0]
        abs_part_error[i] = abs_part_error[i] + abs(data[index,0])
        number_of_points = number_of_points + 1
    print(number_of_points)
    part_error[i] = part_error[i]/number_of_points
    abs_part_error[i] = abs_part_error[i]/number_of_points
    

#abs_part_error[:11] = part_error[:11] 
plt.plot(divs, part_error, linewidth=6, color = "darkorange")
plt.xlabel("Divisions relative to the loss", fontsize=25)
plt.ylabel("Partitioning error", fontsize=25)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

#plt.yticks([-0.15, -0.10, -0.05,  0, 0.05], fontsize = 17)
plt.xlim([-9, 7])
plt.ylim([-0.2, 0.2])
plt.xticks([-5, 0, 5])
plt.yticks([-0.2, -0.1, 0, 0.1, 0.2])

#plt.savefig("fig-part.pdf")
    
