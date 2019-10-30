# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:41:56 2019

@author: Abdul Hannan
"""

import numpy as np
import scipy.linalg as sla
from Random_Orthogonal import OrthVec

class SIM:
    def __init__(self, n,m):
        self.n = n
        self.m = m
        self.a = np.zeros((self.n, self.n))
        for i in range (self.n):
           self.a[i][i] = np.random.randint(1, 100)
        print (self.a)
        self.B= np.eye(self.n)
#        self.a = np.array([[2,np.sqrt(2),4],[np.sqrt(2),6,np.sqrt(2)],[4,np.sqrt(2),2]])
#        self.a = [[0,1,0,0,0,0],[1,0,1,0,0,0],[0,1,0,1,0,0],[0,0,1,0,1,0],[0,0,0,1,0,1],[0,0,0,0,1,0]]
        SV = OrthVec(self.n,self.m)
#        SV = OrthVec(self.n,self.m)
        SV.Calc()
        self.R = SV.vt
        
    def Calculate(self):
        F= open("Eigen.csv",'w')
        eig = np.zeros(self.m)
        err = 1
        while (err> 1e-20):
            Temp = np.linalg.solve(self.a,self.R) #Solve for Th using A*Th = R
            Th = Temp/(np.linalg.norm(Temp))
            Tht = np.transpose(Th)
            Ax = np.matmul(Tht,self.a)  # Reorientation
            Ax = np.matmul(Ax,Th)
#            Bx = np.matmul(Tht,self.B)
#            Bx = np.matmul(Bx,Th)
#            L = sla.qz(Ax, Bx)
#            self.y = (np.diag(L[0])/np.diag(L[1]))
            y = np.linalg.eig(Ax)       #Calculation of eigenvaluesand eigenvectors
            y= list(y)
            self.phi = y.pop(1)
            self.R = np.matmul(Th, self.phi)
            self.y = np.ravel(y)
            F.write(str(self.y))
            F.write("\n")
            err = max(abs(self.y - eig))
            eig = self.y                
        print(self.y)
        print()
        F.close()

if __name__ == '__main__':
    Test = SIM(3,2)
    Test.Calculate()
