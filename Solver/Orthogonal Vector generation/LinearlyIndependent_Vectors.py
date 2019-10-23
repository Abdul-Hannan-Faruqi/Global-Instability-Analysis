# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:38:35 2019

@author: Abdul Hannan
"""

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

size = 3
v = []
v1 = np.random.rand(size)
v1n = v1/np.linalg.norm(v1)
v.append(v1n)

i=1

while (i<size):
    vr = np.random.rand(size)
    vrmag = np.linalg.norm(vr)
    vn = vr/vrmag
    v0 = vr
    v0mag = vrmag
    
    for vector in v:
        v0mag = v0mag - np.dot(v0, vector)/np.linalg.norm(vector)
    if v0mag > 1e-8* vrmag:
        v.append(list(vn))
        i+=1

v= np.array(v)
plv = np.ndarray(shape = (3,6))        
for i in range (len(v)):
    plv[i] = list(np.zeros(size))+ list(v[i])
X,Y,Z,P,Q,R = zip(*plv) 
   
if __name__ == "__main__":
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
#    ax.plot3D(v[:,0], v[:,1], v[:,2], 'k.', alpha = 0.8)
    ax.quiver(X,Y,Z,P,Q,R)
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    plt.show()
