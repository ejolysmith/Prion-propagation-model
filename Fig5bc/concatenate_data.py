# -*- coding: utf-8 -*-
"""
This code concatenates the outputs from prion-model_big-sweep.c 
into a single txt file called data.txt
"""

import numpy as np
import numpy


c1 = 'prion-data-model-6-bigsweep-lam1p75-uni-2-3-v1.txt'
c2 = 'prion-data-model-6-bigsweep-lam1p75-uni-2-3-v2.txt'
c3 = 'prion-data-model-6-bigsweep-lam1p75-uni-3-2-v1.txt'
c4 = 'prion-data-model-6-bigsweep-lam1p75-uni-3-3-v1.txt'
c5 = 'prion-data-model-6-bigsweep-lam1p75-uni-3-3-v2.txt'
c6 = 'prion-data-model-6-bigsweep-lam1p75-uni-3-4-v1.txt'
c7 = 'prion-data-model-6-bigsweep-lam1p75-uni-3-4-v2.txt'
c8 = 'prion-data-model-6-bigsweep-lam1p75-uni-4-1-v1.txt'
c9 = 'prion-data-model-6-bigsweep-lam1p75-uni-4-2-v1.txt'
c10 = 'prion-data-model-6-bigsweep-lam1p75-uni-4-2-v2.txt'

data1 = np.loadtxt(c1)
data2 = np.loadtxt(c2)
data3 = np.loadtxt(c3)
data4 = np.loadtxt(c4)
data5 = np.loadtxt(c5)
data6 = np.loadtxt(c6)
data7 = np.loadtxt(c7)
data8 = np.loadtxt(c8)
data9 = np.loadtxt(c9)
data10 = np.loadtxt(c10)


data = np.concatenate((data1, data2, data3, data4, data5, data6, data7, data8, data9, data10), axis = 0)


        
np.savetxt('data.txt', data)