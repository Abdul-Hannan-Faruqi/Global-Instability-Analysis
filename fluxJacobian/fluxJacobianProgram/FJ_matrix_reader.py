# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 16:54:49 2020

@author: Ahmad Faraz
"""

import numpy as np

filename = 'FJ_matrixT.csv'
f = open(filename, 'r')

n = input("Enter size of matrix: ")
n = int(n)
M = np.zeros((n, n))

s = f.readline() #reads header
s = f.readline()
while(s): #while the line read from the file is not empty
    u = s.split(',') #splitting into a list 
    u[0] = int(float(u[0])) #conversion of i_index
    u[1] = int(float(u[1])) #conversion of j_index
    u[2] = float(u[2].rstrip('\n')) #conversion of value at i,j
    M[u[0]][u[1]] = u[2]
    s = f.readline()
    
f.close()