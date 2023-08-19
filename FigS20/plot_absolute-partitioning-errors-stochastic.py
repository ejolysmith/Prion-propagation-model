"""
Generates figS20d
"""

import numpy as np
from matplotlib import pylab as plt



path = 'data_divs-rel-to-loss/'
c1 = path + "prion-divisions-relative-to-loss-0p1.txt"
c2 = path + "prion-divisions-relative-to-loss-1.txt"
c3 = path + "prion-divisions-relative-to-loss-10.txt"


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

part_error1 = np.zeros(len(divs))
abs_part_error1 = np.zeros(len(divs))

for i in range(len(divs)):
    where_i = np.where(data[:,1] == divs[i])[0]
    number_of_points = 0    
    for j in range(len(where_i)):
        index = int(where_i[j])
        part_error1[i] = part_error1[i] + data[index,0]
        abs_part_error1[i] = abs_part_error1[i] + abs(data[index,0])
        number_of_points = number_of_points + 1
    print(number_of_points)
    part_error1[i] = part_error1[i]/number_of_points
    abs_part_error1[i] = abs_part_error1[i]/number_of_points
    
    

data = np.loadtxt(c2)


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
        
     
divs = [-10, -9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,7,8,9,10]

part_error2 = np.zeros(len(divs))
abs_part_error2 = np.zeros(len(divs))

for i in range(len(divs)):
    where_i = np.where(data[:,1] == divs[i])[0]
    number_of_points = 0    
    for j in range(len(where_i)):
        index = int(where_i[j])
        part_error2[i] = part_error2[i] + data[index,0]
        abs_part_error2[i] = abs_part_error2[i] + abs(data[index,0])
        number_of_points = number_of_points + 1
    print(number_of_points)
    part_error2[i] = part_error2[i]/number_of_points
    abs_part_error2[i] = abs_part_error2[i]/number_of_points
    
    
    
    
data = np.loadtxt(c3)


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
        
     
divs = [-10, -9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,7,8,9,10]

part_error3 = np.zeros(len(divs))
abs_part_error3 = np.zeros(len(divs))

for i in range(len(divs)):
    where_i = np.where(data[:,1] == divs[i])[0]
    number_of_points = 0    
    for j in range(len(where_i)):
        index = int(where_i[j])
        part_error3[i] = part_error3[i] + data[index,0]
        abs_part_error3[i] = abs_part_error3[i] + abs(data[index,0])
        number_of_points = number_of_points + 1
    print(number_of_points)
    part_error3[i] = part_error3[i]/number_of_points
    abs_part_error3[i] = abs_part_error3[i]/number_of_points



#abs_part_error[:11] = part_error[:11] 
plt.plot(divs, abs_part_error1, linewidth=6, color = 'b',label='Large fluctuations')
plt.plot(divs, abs_part_error2, linewidth=6, color = 'r',label='Moderate fluctuations')
plt.plot(divs, abs_part_error3, linewidth=6, color = 'g',label='Small fluctuations')

#plt.plot(divs, part_error3, linewidth=6, color = 'k',label='Average chaperone abundance = 100')
plt.xlabel("Divisions relative to the loss", fontsize=16)
plt.ylabel("Partitioning error", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.legend()
#plt.yticks([-0.15, -0.10, -0.05,  0, 0.05], fontsize = 17)
plt.xlim([-8, 7])
plt.ylim([-0.2, 0.2])
plt.xticks([-5, 0, 5])
#plt.yticks([-0.2, -0.1, 0, 0.1, 0.2])
plt.ylim([0,0.2])
plt.yticks([0, 0.05, 0.1, 0.15, 0.2])
plt.savefig("fig-abs-part-stochastic.pdf")
    
