# -*- coding: utf-8 -*-
"""

@author: Ahmad Faraz
"""

"""

Flux Jacobian 005 for natural convection 
(including boundary conditions)

"""

import os
import numpy as np
import pyvista as pv
from findNeighbours import findNeighbours
import naturalConvection_dash_equations005 as eq
import scaleCalc
import dxCalc

"""
Usage:
    Flux Jacobian for Natural Convection
    Change the names of the vtk files
    Change Ra
"""
Ra_string = '10000'
time_string = '50000'
simulationType = 'Simple'

baseFolder = 'test/OnsetCavity/50x50/Ra'+Ra_string+'-'+simulationType+'/VTK/'
gridExtension = 'buoyant'+simulationType+'Foam-uniformMesh-p_rghBC-50x50-'
gridFileName = baseFolder+gridExtension+'Ra'+ Ra_string + '_'+ time_string + '.vtk'
boundaryFileName_1 = baseFolder+'cavityBottomWall/cavityBottomWall_'+time_string+'.vtk'
boundaryFileName_2 = baseFolder+'cavityWalls/cavityWalls_' +time_string+ '.vtk'
boundaryFileName_3 = baseFolder+'atmosphereSurface/atmosphereSurface_'+time_string+'.vtk'
boundaryFileName_4 = baseFolder+'freeAtmosphere/freeAtmosphere_'+time_string+'.vtk'

Ra = 10000
Pr = 0.7116
rho = 1.19

Re1 = (Ra/Pr)**(1/2)
Re2 = (Ra*Pr)**(1/2)

#To non-dimensionalisse results:
U_scale = scaleCalc.U_scale(Ra)
L_scale = scaleCalc.L_scale(Ra)
p_scale = scaleCalc.p_scale(Ra)

print("Reading VTK files ...")
grid = pv.read(gridFileName)
boundaryGrid_1 = pv.read(boundaryFileName_1)
boundaryGrid_2 = pv.read(boundaryFileName_2)
boundaryGrid_3 = pv.read(boundaryFileName_3)
boundaryGrid_4 = pv.read(boundaryFileName_4)
boundaryGrid = boundaryGrid_1 + boundaryGrid_2 + boundaryGrid_3 + boundaryGrid_4
boundaryList = [boundaryGrid_1['patchID'][0],
                boundaryGrid_2['patchID'][1],
                boundaryGrid_3['patchID'][2],
                boundaryGrid_4['patchID'][3]]

cellID = grid.cell_arrays['cellID']
#Contains only the cellIDs of internal cells

n = len(cellID)
#n = Number of cells

boundaryCheckIndex = len(cellID)#All boundaries, must be changed for selected boundaries

points = np.concatenate([grid.cell_centers().points, boundaryGrid.cell_centers().points])
#Contains centers of cells as well as boundary patches
points = points/L_scale

p = np.concatenate([grid.cell_arrays['p'], boundaryGrid.cell_arrays['p']])
#Contains p values of cells as well as boundary patches
p = p/p_scale

U = np.concatenate([grid.cell_arrays['U'], boundaryGrid.cell_arrays['U']])
#Contains U values of cells as well as boundary patches
U = U/U_scale

theta = np.concatenate([grid.cell_arrays['T'], boundaryGrid.cell_arrays['T']])
#Contains theta values of cells as well as boundary patches
theta = (theta - 300)/20

print("Finding cell neighbours ...")
neighbours = findNeighbours(cellID, points)
neighbours = np.int32(neighbours)
#Contains indices of west, east, south and north neighbours of a cell

dx = dxCalc.dxCalc(points, neighbours)
# Calculates (dx[w], dx[e], dx[w], dx[s], dx[n]) for an interior point

boundaryCheck = dxCalc.boundaryChecker(neighbours, n)
# Returns whether a point has w,e,s,n boundaries or not (p in case of no boundary)

boundaryPatchCheck = dxCalc.boundaryPatchChecker(neighbours, n, boundaryGrid.cell_arrays['patchID'])
# Returns which patch a point has w,e,s,n boundaries to or -1 in case of no boundary

n_var = 4
rowLength = n_var*n

print("Writing FJ_matrix ...")
filename = 'FJ_matrix005.csv'
os.remove(filename)
file = open(filename, 'ab')
np.savetxt(file, np.vstack(('i_index', 'j_index', 'value')).T, fmt = '%s', delimiter = ',')
for i in cellID[0:2500]:
    row_0 = eq.p_dash_equation(i, neighbours[i], p, U, theta, dx, boundaryCheck, boundaryPatchCheck, rowLength, n_var, Re1, Re2)
    row_1 = eq.Ux_dash_equation(i, neighbours[i], p, U, theta, dx, boundaryCheck, boundaryPatchCheck, rowLength, n_var, Re1, Re2)
    row_2 = eq.Uy_dash_equation(i, neighbours[i], p, U, theta, dx, boundaryCheck, boundaryPatchCheck, rowLength, n_var, Re1, Re2)
    row_3 = eq.theta_dash_equation(i, neighbours[i], p, U, theta, dx, boundaryCheck, boundaryPatchCheck, rowLength, n_var, Re1, Re2)
    
    row_0 = row_0[0:10000]
    row_1 = row_1[0:10000]
    row_2 = row_2[0:10000]
    row_3 = row_3[0:10000]
    
    row_set = np.vstack((row_0.T, row_1.T, row_2.T, row_3.T))
    i_index, j_index = np.where(row_set != 0)
    value = row_set[i_index, j_index]
    i_index = i*n_var + i_index
    
    np.savetxt(file, np.vstack((i_index, j_index, value)).T, fmt = '%s', delimiter = ',')
    if i%500 == 0:
        print('500 cells entered:.\n')
file.close()    