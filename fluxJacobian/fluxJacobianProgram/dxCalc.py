# -*- coding: utf-8 -*-
"""
@author: Ahmad Faraz
"""

"""
dxCalc

Calculates (dx[w], dx[e], dx[w], dx[s], dx[n]) for an interior point
"""

import numpy as np

def dxCalc(points, neighbours):
    n = len(neighbours)
    dx = np.zeros(points[neighbours].shape)
    
    for i in np.arange(n):
        dx[i] = points[i] - points[neighbours[i]]
        
    return np.abs(dx)

def boundaryChecker(neighbours, n):
    Bchecker = n*np.ones(neighbours.shape)
    boundarys = neighbours>=Bchecker # Truth values
    boundarysAssigned = boundarys.tolist()
    wesn = ['w', 'e', 's', 'n']
    
    for i in np.arange(n):
        for j in [0, 1, 2, 3]:
            if boundarys[i][j] == True:
                boundarysAssigned[i][j] = wesn[j]
            else:
                boundarysAssigned[i][j] = 'p'     
    return boundarysAssigned

def boundaryPatchChecker(neighbours, n, bPatches):
    boundarysAssigned = np.zeros(neighbours.shape)
    
    for i in np.arange(n):
        for j in [0, 1, 2, 3]:
            if neighbours[i][j]>=n:
                boundarysAssigned[i][j] = bPatches[neighbours[i][j]- n]
            else:
                boundarysAssigned[i][j] = -1
    
    return boundarysAssigned