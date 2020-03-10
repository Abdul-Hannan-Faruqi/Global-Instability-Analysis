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
import time

class Mat_Gen:
    def __init__(self, n, t):
        self.n = n
        self.a = np.zeros((self.n, self.n))        

    def Symmetric(self, max):
        for i in range(self.n):
            j=i
            while j<self.n:
                self.a[i][j] = np.random.randint(1, max)
                if i!= j:
                    self.a[j][i] = self.a[i][j]
                j = j+1
        return self.a
                    
    def Diagonal(self, max):
        for i in range (self.n):
            self.a[i][i] = np.random.randint(1, max)
        return self.a
            
    def Random(self, max):
        for i in range (self.n):
            for j in range(self.n):
                self.a[i][j] = np.random.randint(1, max)#+ 1j*(np.random.randint(1, max))
        return self.a
    
class SIM:
    def __init__(self, n, m):
        self.n = n
        self.m = m
        dsum = 0
        self.B= np.eye(self.n)
        a = np.load("Matrices\M.npy")
        m = np.load("Matrices\A.npy")
#        Matrix = Mat_Gen(self.n, t)
#        a = np.array([[1,2,3],[0,4,5],[0,0,1]])
#        a= 100*np.random.rand(self.n, self.n)
#        a = np.array([[1,4,3,5],[4,3,5,6],[3,5,4,8],[5,6,8,9]])
#        a=[[1,-1],[1,1]]
#        a = np.array([[1,2,3],[2,3,4],[3,4,5]])
#        a=np.array([[0,2,3,4,5,6],[2,0,4,5,6,7],[3,4,0,6,7,8],[4,5,6,0,8,9],[5,6,7,8,0,10],[6,7,8,9,10,11]])
#        a = np.array([[4/5,-3/5,0],[3/5,4/5,0],[1,2,2]])
#        a = np.array([[2,np.sqrt(2),4],[np.sqrt(2),6,np.sqrt(2)],[4,np.sqrt(2),2]])
#        a = np.array([[1,2,2],[0,2,1],[-1,2,2]], dtype = 'complex')
        print (a)
        k = np.linalg.inv(a)@m
        self.k = k
        np.save("Matrices\K.npy", self.k)
        SV = OrthVec(self.n,self.m)
        SV.Calc()
        R = list(SV.vt)
        self.R = np.array(R, dtype= 'complex')
        self.Calculate(k)
        
    def Calculate(self, k):
        F= open("Solution\Eigen.txt",'w')
        err = 1
        eig = 0
        iter = 0
        anew = np.zeros((self.m, self.m))
        start_time = time.process_time()
        while (err > 1e-4):
#            Th = np.linalg.solve(k,self.R) #Solve for Th using A*Th = R  (for smallest eigenvalues)
            Th = k@self.R  #for largest eigenvalues
            Th = Th/ max(np.ravel(Th))
            Th, r = np.linalg.qr(Th)
            Ax = Th.T@k@Th
            y, self.phi = np.linalg.eig(Ax)       # Calculation of eigenvalues and eigenvectors
#            if iter%1 == 0:
#                self.Plot(iter)               
            self.R = np.matmul(Th, self.phi)    #Eigenvector
            self.y = np.ravel(y)                #Eigenvalues
            F.write(str(self.y))
            F.write("\n")            
            err1 = max(abs(k@self.R[:,1])-(self.y[1]*self.R[:,1]))
            err2 = max(abs(eig-self.y))
            err = min(err1, err2)
            if (iter == 100):
                self.Complex_Shift()
            eig = self.y
            iter += 1 
        print("Execution Time = ", (time.process_time()-start_time), "s")   
        self.y = self.y + self.shift         
        print(self.y)
        print("shift = ", self.shift)
        np.save("Solution\R.npy", self.R)
        print("\nIter = ", iter)
        F.close()
        
    def Multigrid(self, Ax,y):
        b = np.ndarray(shape=(self.m, self.m), dtype= 'complex')
        for i in range(self.m):
            for j in range(self.m):
                if i!=j:
                    b[i][j]= Ax[i][j]
                else:
                    b[i][j] = Ax[i][j]-self.y[0]
        print("b=", b)

    def Complex_Shift(self):
        self.shift = min(self.y)
        self.k = self.k - self.shift*(np.eye(self.n))
                   
    def Plot(self, iter):
        plv = np.ndarray(shape = (3,6)) 
        for i in range (3):
            if self.m==2:
                if i==self.m:
                    plv[i] = list(np.zeros(3))+ list(np.zeros(3))
                else:
                    plv[i] = list(np.zeros(3))+ list(self.phi[i])+list(np.zeros(1))
            if self.m==3:
                plv[i] = list(np.zeros(3))+ list(self.R[:,i])
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        X,Y,Z,P,Q,R = zip(*plv) 
        ax.quiver(X,Y,Z,P,Q,R)
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])

        plt.savefig("v"+str(iter)+".jpg")
        plt.close()

if __name__ == '__main__':
    n = input("Enter the dimension of matrix (n): ")
    m = input("Enter the dimension of subspace (m): ")
    Test = SIM(int(n), int(m))    # (dimension of matrix, dimension of sub-space)