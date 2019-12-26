# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:41:56 2019

@author: Abdul Hannan
"""

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sla
from Random_Orthogonal import OrthVec

class SIM:
    def __init__(self, n,m):
        self.n = n
        self.m = m
        self.a = np.zeros((self.n, self.n))
#        for i in range (self.n):
#            self.a[i][i] = np.random.randint(1, 1000)
#            for j in range(self.n):
#                self.a[i][j] = np.random.randint(1, 10)#+ 1j*(np.random.randint(1, 100))
        self.B= np.eye(self.n)
#        self.a= 100*np.random.rand(self.n, self.n)
#        self.a=[[1,-1],[1,1]]
#        self.a = np.array([[1,2,3],[2,3,4],[3,4,5]])
#        self.a=np.array([[0,2,3,4,5,6],[2,0,4,5,6,7],[3,4,0,6,7,8],[4,5,6,0,8,9],[5,6,7,8,0,10],[6,7,8,9,10,11]])
#        self.a = np.array([[4/5,-3/5,0],[3/5,4/5,0],[1,2,2]])
#        self.a = np.array([[2,np.sqrt(2),4],[np.sqrt(2),6,np.sqrt(2)],[4,np.sqrt(2),2]])
#        self.a = np.array([[1,2,2],[0,2,1],[-1,2,2]])
        print (self.a)
        SV = OrthVec(self.n,self.m)
        SV.Calc()
        self.R = SV.vt
        
    def Calculate(self):
        F= open("Eigen.txt",'w')
        eig = np.zeros(self.m)
        err = 1
        iter = 0
        self.a= self.a/np.linalg.norm(self.a)
        while (err> 1e-4):
            Th = np.linalg.solve(self.a,self.R) #Solve for Th using A*Th = R
#            Th = np.matmul(self.a, self.R)
            Th = Th/(np.linalg.norm(Th))
            if np.dot(Th[:,0], Th[:,1])< 1e-5:
               Th, r = np.linalg.qr(Th)
            Tht = np.transpose(Th)
            Ax = np.matmul(Tht,self.a)  # Reorientation
            Ax = np.matmul(Ax,Th)
#            Bx = np.matmul(Tht,self.B)
#            Bx = np.matmul(Bx,Th)
            y = np.linalg.eig(Ax)       #Calculation of eigenvalues and eigenvectors
            y= list(y)
            self.phi = y.pop(1)
#            self.Plot(iter)               
            self.R = np.matmul(Th, self.phi)    #Eigenvector
            self.y = np.ravel(y)                #Eigenvalues
            for i in range(self.m):
                F.write(str(self.y))
                F.write(',')
            F.write("\n")
            err = max(abs(self.y - eig))
            eig = self.y
            iter+= 1               
        print(self.y)
        print()
        F.close()
         
    def Plot(self, iter):
        plv = np.ndarray(shape = (3,6)) 
        for i in range (3):
            if self.m==2:
                if i==self.m:
                    plv[i] = list(np.zeros(3))+ list(np.zeros(3))
                else:
                    plv[i] = list(np.zeros(3))+ list(self.phi[i])+list(np.zeros(1))
            if self.m==3:
                plv[i] = list(np.zeros(3))+ list(self.R[:,i])#self.phi[:,i])
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        X,Y,Z,P,Q,R = zip(*plv) 
        ax.quiver(X,Y,Z,P,Q,R)
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])

        plt.savefig("b"+str(iter)+".jpg")
        plt.close()

if __name__ == '__main__':
    Test = SIM(3,3)
    Test.Calculate()
