# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 20:14:30 2022

@author: Iván
"""

import numpy as np
import matplotlib.pyplot as plt

n_1=np.linspace(0, 0.28, 100)

n_2=np.linspace(0.28, 1, 100)

T=1

gamma=10

lb2=1/gamma

P_1=T*n_1

P_2=T*(n_2*gamma*lb2*n_2**2/3)
        
plt.plot(n_1, P_1, color='blue')
plt.plot(n_2, P_2, color='blue')
plt.vlines(0.28, P_2[0], P_1[-1], ls='dashed')
plt.grid(True)
plt.xlabel(r'$n$ [rods/m$^3$]')
plt.ylabel(r'$P$ [Pa]')
plt.title(r'Pressure ($P$) as a function of $n$')
plt.savefig('P_vs_n.pdf')
