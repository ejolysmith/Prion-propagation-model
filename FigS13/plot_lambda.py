#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code creates Fig S13. 
The C program compiled from prion-model_lambda is run for different lambda values, 
which are set manually in the code on line 53. The output is an array of partitioning errors, 
we take the average absolute partitioning errors and save them here for each lambda along with 
their error taken as the standard error of the mean. 
"""
import numpy as np
from matplotlib import pylab as plt
from scipy.optimize import curve_fit

measured_part = 0.051
lamb = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,  1, 1.5, 1.6, 1.7, 1.8, 1.9])*50/np.log(2) #Avg protein concentration
part = [0.229, 0.153, 0.124, 0.10745, 0.0944, 0.088, 0.0790, 0.0762, 0.0705, 0.0660, 0.0541, 0.0532, 0.0516, 0.0504, 0.0483] #Avg abs part error
err = [0.002, 0.001, 0.0007, 0.0008, 0.0009, 0.001, 0.0005, 0.0006,  0.0007, 0.0005, 0.0003, 0.0004, 0.0004, 0.0003, 0.0006] #Error on avg abs part error


def func(x, a):
    return a/np.sqrt(x)
popt, pcov = curve_fit(func, lamb, part, sigma = err)

x = np.linspace(0, max(lamb), 100)
plt.plot(x, func(x, *popt), '-', label = r'Fitted function of the form y = $A/\sqrt{x}$', color = 'C1')
plt.plot(x, measured_part*np.ones(len(x)), '--', label = "Measured partitioning error", color = 'C0')
plt.plot(lamb, part,  '.', label = 'Simulations', color = 'k', markersize = 8)
plt.plot([126.2], func(126.2,*popt), '.', color = 'r', label = 'Point of intersect, $\lambda = 1.75$', markersize = 15)
plt.xlim([min(lamb)-3, max(lamb)+3])
plt.ylim([0,0.3])
plt.xlabel(r'Average protein concentration $\lambda \tau_{div} / \ln(2)$')
plt.ylabel('Average absolute partitioning error')
plt.legend()

plt.savefig('fig-SI-lambda.pdf')