# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:38:35 2019

@author: Abdul Hannan
"""

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

class OrthVec():
    def __init__(self, size, num):
        self.size = size
        self.num = num
        self.v = []
        v1 = np.random.rand(size)
        v1n = v1/np.linalg.norm(v1)
        self.v.append(list(v1n))
        
    def Calc(self):
        i=1
        while (i<self.num):
            v0 = np.random.rand(self.size)    
            for vector in self.v:
                v0 = v0 - ((np.dot(v0, vector)/np.dot(vector, vector))*np.array(vector))    
            vn = v0/np.linalg.norm(v0)
            self.v.append(list(vn))
            i+=1
            
        self.v = np.array(self.v)
        self.vt = np.transpose(self.v)

   
if __name__ == "__main__":
    Xa = OrthVec(6,3)
    Xa.Calc()
    vecs= Xa.v
    plv = np.ndarray(shape = (Xa.size,2*Xa.size))        
    for i in range (Xa.size):
        plv[i] = list(np.zeros(Xa.size))+ list(vecs[i])
    X,Y,Z,P,Q,R = zip(*plv) 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
#    ax.plot3D(v[:,0], v[:,1], v[:,2], 'k.', alpha = 0.8)
    ax.quiver(X,Y,Z,P,Q,R)
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    plt.show()
