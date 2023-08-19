#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: euan

This code takes in the simulated time-of-loss data from 
the compiled prion-model_big-sweep.c, which simulates the 
stochastic model multiple times with different parameters 
and saves the time-of-loss at each point. Here this code 
turns the scattered time-of-loss measurements into a 
contour plot.  

"""



import numpy as np
from matplotlib import pylab as plt
import numpy

data = np.loadtxt('data.txt')


x = []
y = []
z = []
for i in range(len(data[:,0])):
    if data[:,7][i] > 0 and data[:,6][i] > 0:
        x.append((data[:,7][i]))
        y.append((data[:,6][i]))
        z.append((data[:,0][i]) )
        
        
  
    
  
 #These two points are added to the data to make the top right region dark blue
 #There are no actual data points from the top right region due to long simulations 
 #times as indicated in Fig. S14. 
  
x.append(0.9*1e-2)
y.append(0.9*1e-2)
z.append(5000)


x.append(0.9*1e-1)
y.append(0.9*1e-1)
z.append(5000)






x = np.array(x)
y = np.array(y)
z = np.array(z)








xy = np.vstack([x,y])


#Binning data and taking averages to smooth the contours
num_bins = 16
x_bined_exponent = np.linspace(-5, -1, num_bins)
y_bined_exponent = np.linspace(-5, -1, num_bins)
x_bined = np.zeros(len(x_bined_exponent))
y_bined = np.zeros(len(y_bined_exponent))
z_bined = np.zeros((len(y_bined), len(y_bined)))
z_bined_numbers = np.zeros((len(y_bined), len(y_bined)))

for i in range(num_bins):
    x_bined[i] = 10**x_bined_exponent[i]
    y_bined[i] = 10**y_bined_exponent[i]


total_numbers = 0

x_bined_3 = np.zeros((len(y_bined), len(y_bined)))
y_bined_3 = np.zeros((len(y_bined), len(y_bined)))

for i in range(len(x)):
    a = 0
    for j in range(num_bins - 1):
        if x_bined[j] <= x[i] < x_bined[j+1]:
            for k in range(num_bins-1):
                if y_bined[k] <= y[i] < y_bined[k+1] :
                    z_bined[k,j] = z_bined[k,j] + z[i]
                    z_bined_numbers[k,j] = z_bined_numbers[k,j] + 1
                    a = 1
                    total_numbers = total_numbers + 1
                    break
            #j = j + 1
for i in range(len(x_bined)-1):
    x_bined[i] = x_bined[i+1]*0.5 + x_bined[i]*0.5
    y_bined[i] = y_bined[i+1]*0.5 + y_bined[i]*0.5

        
for j in range(num_bins):
    for k in range(num_bins):
        if z_bined_numbers[j,k] > 0:
            z_bined[j,k] = z_bined[j,k]/z_bined_numbers[j,k]
            
x_bined_2 = []
y_bined_2 = []
z_bined_2 = []
for j in range(num_bins):
    for k in range(num_bins):
        if z_bined_numbers[k,j] > 0 :
            x_bined_2.append(x_bined[j])
            y_bined_2.append(y_bined[k])
            z_bined_2.append(z_bined[k,j])
 




#Plotting
fig, (ax1) = plt.subplots(nrows=1)



levels = [0,100,200,300,400,500,1000]
density = plt.tricontourf(x_bined_2,y_bined_2,z_bined_2, levels = levels, cmap="Blues", vmin = 0, vmax = 1000, extend='max')#"RdBu_r")
plt.tricontour(x_bined_2,y_bined_2,z_bined_2, alpha = 0.5, levels = levels, linewidths=0.5, colors='k', vmin = 0, vmax = 1000)
CS = plt.tricontour(x_bined_2,y_bined_2,z_bined_2, levels = [129.26], linewidths=2, colors='darkorange', vmin = 0, vmax = 1000)

alphadot = 1.97000000e-04
gammadot = 1.12200000e-03
plt.scatter([gammadot], [alphadot], color = 'darkorange', s = 200)
plt.xscale('log')
plt.yscale('log')
cb = plt.colorbar(density, pad=0.01)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)

cb.set_label(label='Average time of loss (min)',size = 15, rotation=-90, labelpad = 17)
cb.set_ticklabels([0, 100, 200, 300, 400, 500, '>1000'])

plt.xlabel(r"Elongation rate $\gamma V_{0}$ (min$^{-1}$)", fontsize = 15)
plt.ylabel(r"Fragmentation rate $\alpha V_{0}$ (min$^{-1}$)", fontsize = 15)

plt.plot( np.linspace(gammadot, 2e-3, 10), np.linspace(alphadot, 1e-3, 10), color = 'darkorange', linewidth = 3)
plt.scatter(np.linspace(gammadot, 2e-3, 10)[4], np.linspace(alphadot, 1e-3, 10)[4], color = 'darkorange', s = 50)
plt.scatter(np.linspace(gammadot, 2e-3, 10)[6], np.linspace(alphadot, 1e-3, 10)[6], color = 'darkorange', s = 50)
plt.scatter(np.linspace(gammadot, 2e-3, 10)[7], np.linspace(alphadot, 1e-3, 10)[7], color = 'darkorange', s = 50)
plt.scatter(np.linspace(gammadot, 2e-3, 10)[8], np.linspace(alphadot, 1e-3, 10)[8], color = 'darkorange', s = 50)


#plt.savefig('figure-time-contour.pdf')
