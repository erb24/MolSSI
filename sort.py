#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 20:52:21 2017

@author: ebeyerle

sort.py--Code for sorting LE4PD modes by relaxation time.
"""

import os
import numpy as np

pwd=os.getcwd()
tau_indx=np.loadtxt(pwd+"/tau_scaled.dat")
tau=tau_indx[:,1]
n=len(tau)
tau_sorted=np.sort(tau)
tau_sorted_new=np.zeros(n)
for i in range(0,n):
    tau_sorted_new[i]=tau_sorted[n-1-i]

dummy=np.zeros((n,2))
for i in range(0,n): 
    c=0   
    while (c <=n-1):
        if tau_sorted_new[i] == tau_indx[c,1]:
            dummy[i,1],dummy[i,0]=tau_sorted_new[i],tau_indx[c,0]
        c=c+1
        
np.savetxt(pwd+'/tau_sorted.dat',dummy)
