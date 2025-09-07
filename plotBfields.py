#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 13:35:15 2020

@author: bret
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit



data = np.loadtxt(fname = 'Bfields.txt', delimiter =None, skiprows = 0)

plt.plot(data[:,0],data[:,1])
plt.xlabel('z (cm)')
plt.ylabel('Bx (G)')
plt.show()

plt.plot(data[:,0],data[:,2])
plt.xlabel('z (cm)')
plt.ylabel('By (G)')
plt.show()

plt.plot(data[:,0],data[:,3])
plt.xlabel('z (cm)')
plt.ylabel('Bz (G)')
plt.show()


