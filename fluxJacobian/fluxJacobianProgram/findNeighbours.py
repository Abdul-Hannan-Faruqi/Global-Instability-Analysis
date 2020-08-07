# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 10:30:33 2020

@author: Ahmad Faraz
"""

import numpy as np
import scipy.spatial as spsp

#Function to find the indices of neighbours
def findIndex(cellID, kNeighbours, points, direction):
    
    if direction == 'w':
        for i in kNeighbours:
            if points[cellID][0] - points[i][0] > 0:
                return i
            
    if direction == 'e':
        for i in kNeighbours:
            if points[cellID][0] - points[i][0] < 0:
                return i
            
    if direction == 's':
        for i in kNeighbours:
            if points[cellID][1] - points[i][1] > 0:
                return i
            
    if direction == 'n':
        for i in kNeighbours:
            if points[cellID][1] - points[i][1] < 0:
                return i
            
    return 'increase k'

#Function to find the neighbours of all cells
def findNeighbours(cellID, points):
    pointsTree = spsp.KDTree(points)
    #Creates a spatial tree for all the points

    neighbours = np.zeros((len(cellID), 4))
    #Initialises a matrix that would contain the east, west, south and north indices of a cell
    
    for i in cellID:
        k = 5
        flag = 0
        while(flag == 0):
            kNeighbours = pointsTree.query(points[i], k)[1]
            #Returns indices of k=5 closest neighbours of the point at index i
        
            indexPosition = -1
            for j in ['w', 'e', 's', 'n']:
                index = findIndex(i, kNeighbours, points, j)
                if  index == 'increase k':
                    k = k + 1
                    break
                else:
                    indexPosition = indexPosition + 1
                    neighbours[i][indexPosition] = index
            
            if indexPosition == 3:
                flag = 1
        print('Cell number: ', i, ', Number of k iterations = ', k-5)
    return neighbours # (w, e, s, n) format
                