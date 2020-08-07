# -*- coding: utf-8 -*-
"""
@author: Ahmad Faraz
"""

"""
Natural Convection Equations 005
"""

import numpy as np

def dxCall(dx):
 #Assignment of various dx
    dx_w = dx[0][0][0]
    dx_e = dx[1][0][0]
    dx_s = dx[2][1][0]
    dx_n = dx[3][1][0]
    
    deltaX = dx[0][0][0] + dx[1][0][0]
    deltaY = dx[2][1][0] + dx[3][1][0]

    return dx_w, dx_e, dx_s, dx_n, deltaX, deltaY

def pWestFunction(a,b,c,d,wPatch):
    if wPatch==3 or wPatch==2 or wPatch==1:
        b[0] = 0
    if wPatch==0:
        b[0] = 0
    return a,b,c,d

def pEastFunction(a,b,c,d,ePatch):
    if ePatch==3 or ePatch==2 or ePatch==1:
        b[1] = 0
    if ePatch==0:
        b[1] = 0
    return a,b,c,d

def pSouthFunction(a,b,c,d,sPatch):
    if sPatch==3 or sPatch==2 or sPatch==1:
        c[2] = 0
    if sPatch==0:
        c[2] = 0
    return a,b,c,d

def pNorthFunction(a,b,c,d,nPatch):
    if nPatch==3 or nPatch==2 or nPatch==1:
        c[3] = 0
    if nPatch==0:
        c[3] = 0
    return a,b,c,d

def p_dash_equation(index,neighbours,p,U,theta,dx,wesn,bPatches,n,n_var,Re1,Re2):
    
    row = np.zeros((n_var*n, 1))
    directions = 4
    a = np.zeros((directions+1,1))#Coefficients:p in wesnp order
    b = np.zeros((directions+1,1))#Coefficients:Ux
    c = np.zeros((directions+1,1))#Coefficients:Uy
    d = np.zeros((directions+1,1))#Coefficients:theta
    
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    
    b[0] = 1/(dx_w + dx_e)
    b[1] = -b[0]
    c[2] = 1/(dx_s + dx_n)
    c[3] = -c[2]
    
    #Assigning coefficients considering boundary points
    if wesn[0]=='w':
        a,b,c,d = pWestFunction(a,b,c,d,Re1,dx,bPatches[0])
    if wesn[1]=='e':
        a,b,c,d = pEastFunction(a,b,c,d,Re1,dx,bPatches[1])
    if wesn[2]=='s':
        a,b,c,d = pSouthFunction(a,b,c,d,bPatches[2])
    if wesn[3]=='n':
        a,b,c,d = pNorthFunction(a,b,c,d,bPatches[3])
    
    k=0
    for i in neighbours:#Among all the neighbours,
        if i<n: #if the neighbour is a inerior point
            row[i*n_var + 1] = b[k]
            row[i*n_var + 2] = c[k]
        k = k + 1
    
    return row

def UxWestFunction(a,b,c,d,Re1,dx,wPatch):
    #Function to assigne coefficients to points wich are adjacent to a western 
    #boundary for different boundary types
    A = np.zeros(a.shape)
    B = np.zeros(b.shape)
    
    #Assignment of varios dx
    dx_w = dx[0][0]
    dx_e = dx[1][0]
    dx_s = dx[2][1]
    dx_n = dx[3][1]
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    #Assignment of coefficients in the pressure boundary condition at no slip
    #walls
    A[1] = -dx_w/(dx_w + dx_e)
    A[4] = (dx_w + dx_e)/dx_w
    A[0] = -(A[1] + A[4])
    B[1] = -2/(Re1*(dx_w + dx_e))
    B[4] = 2/(Re1*(dx_w))
    
    #Assignement of coefficients on no slip walls
    if wPatch==3 or wPatch==2 or wPatch==1 or wPatch==0:
        a[4] = -(a[0]/A[0])*A[4]
        a[1] = a[1] - (a[0]/A[0])*A[1]
        b[0] = 0
        b[1] = b[1] - (a[0]/A[0])*B[1]
        b[4] = b[4] - (a[0]/A[0])*B[4]
        a[0] = 0
        
    return a,b,c,d

def UxEastFunction(a,b,c,d,Re1,dx,ePatch):
    #Function to assigne coefficients to points wich are adjacent to an eastern 
    #boundary for different boundary types
    A = np.zeros(a.shape)
    B = np.zeros(b.shape)
    
    #Assignment of varios dx
    dx_w = dx[0][0]
    dx_e = dx[1][0]
    dx_s = dx[2][1]
    dx_n = dx[3][1]
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    #Assignment of coefficients in the pressure boundary condition at no slip
    #walls
    A[0] = dx_e/(dx_w + dx_e)
    A[4] = -(dx_w + dx_e)/dx_e
    A[1] = -(A[0] + A[4])
    B[0] = -2/(Re1*(dx_w + dx_e))
    B[4] = 2/(Re1*(dx_e))
    
    #Assignement of coefficients on no slip walls
    if ePatch==3 or ePatch==2 or ePatch==1 or ePatch==0:
        a[4] = -(a[1]/A[1])*A[4]
        a[0] = a[0] - (a[1]/A[1])*A[0]
        b[1] = 1
        b[0] = b[0] - (a[1]/A[1])*B[0]
        b[4] = b[4] - (a[1]/A[1])*B[4]
        a[1] = 1
        
    return a,b,c,d

def UxSouthFunction(a,b,c,d,dx,sPatch):
    if sPatch==3 or sPatch==2 or sPatch==1:
        b[2] = 0
    if sPatch==0:
        
        #Assignment of varios dx
        dx_w = dx[0][0]
        dx_e = dx[1][0]
        dx_s = dx[2][1]
        dx_n = dx[3][1]
        dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
        phi_N = -(dx_s)/(dx_n*(dx_s+dx_n))
        phi_P = (dx_s+dx_n)/(dx_n*dx_s)
        phi_S = -(phi_N + phi_P)
        
        b[3] = b[3] - b[2]*(phi_N/phi_S)
        b[4] = b[4] - b[2]*(phi_P/phi_S)
        b[2] = 0
    return a,b,c,d

def UxNorthFunction(a,b,c,d,dx,nPatch):
    if nPatch==3 or nPatch==2 or nPatch==1:
        b[3] = 0
    if nPatch==0:
        
        #Assignment of varios dx
        dx_w = dx[0][0]
        dx_e = dx[1][0]
        dx_s = dx[2][1]
        dx_n = dx[3][1]
        dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
        phi_S = (dx_n)/(dx_s*(dx_s+dx_n))
        phi_P = -(dx_s+dx_n)/(dx_n*dx_s)
        phi_N = -(phi_S + phi_P)
        
        b[2] = b[2] - b[3]*(phi_S/phi_N)
        b[4] = b[4] - b[3]*(phi_P/phi_N)
        b[3] = 0
    return a,b,c,d

def Ux_dash_coefficients(U,dx,Re1,neighbours,index,wesn,bPatches):
    directions = 4
    a = np.zeros((directions+1,1))#Coefficients:p in wesnp order
    b = np.zeros((directions+1,1))#Coefficients:Ux
    c = np.zeros((directions+1,1))#Coefficients:Uy
    d = np.zeros((directions+1,1))#Coefficients:theta
    
    #deltaX = dx[0][0]+dx[1][0]
    #deltaY = dx[2][1]+dx[3][1]
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    #Assigning coefficients for interior points without considering boundary
    #points
    a[0] = 1/deltaX
    b[0] = 2/(Re1*dx_w*deltaX) + U[index][0]/deltaX
    a[1] = -1/deltaX
    b[1] = 2/(Re1*dx_e*deltaX) - U[index][0]/deltaX
    b[2] = 2/(Re1*dx_s*deltaY) + U[index][1]/deltaY
    b[3] = 2/(Re1*dx_n*deltaY) - U[index][1]/deltaY
    b[4] = -2/(Re1*dx_e*dx_w) - 2/(Re1*dx_s*dx_n) - (U[neighbours[1]][0]-U[neighbours[0]][0])/deltaX
    c[4] = - (U[neighbours[3]][0]-U[neighbours[2]][0])/deltaY
    
    #Assigning coefficients considering boundary points
    if wesn[0]=='w':
        a,b,c,d = UxWestFunction(a,b,c,d,Re1,dx,bPatches[0])
    if wesn[1]=='e':
        a,b,c,d = UxEastFunction(a,b,c,d,Re1,dx,bPatches[1])
    if wesn[2]=='s':
        a,b,c,d = UxSouthFunction(a,b,c,d,dx,bPatches[2])
    if wesn[3]=='n':
        a,b,c,d = UxNorthFunction(a,b,c,d,dx,bPatches[3])

    return a,b,c,d


def Ux_dash_equation(index,neighbours,p,U,theta,dx,wesn,bPatches,n,n_var,Re1,Re2):
    
    row = np.zeros((n_var*n, 1))

    a,b,c,d = Ux_dash_coefficients(U,dx,Re1,neighbours,index,wesn,bPatches)
    
    k=0
    for i in np.append(neighbours, index):#Among all the neighbours,
        if i<n: #if the neighbour is a inerior point
            row[i*n_var + 0] = a[k]
            row[i*n_var + 1] = b[k]
            row[i*n_var + 2] = c[k]
            row[i*n_var + 3] = d[k]
        k = k + 1
    
    return row

def UyWestFunction(a,b,c,d,Re1,dx,wPatch):
    
    #Assignement of coefficients on no slip walls
    if wPatch==3 or wPatch==2 or wPatch==1:
        c[0] = 0
    if wPatch==0:
        
        #Assignment of varios dx
        dx_w = dx[0][0]
        dx_e = dx[1][0]
        dx_s = dx[2][1]
        dx_n = dx[3][1]
        dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
        phi_E = -(dx_w)/(dx_e*(dx_w + dx_e))
        phi_P = (dx_e + dx_w)/(dx_e*dx_w)
        phi_W = -(phi_E + phi_P)
        
        c[1] = c[1] - c[0]*(phi_E/phi_W)
        c[4] = c[4] - c[0]*(phi_P/phi_W)
        c[0] = 0    
    return a,b,c,d

def UyEastFunction(a,b,c,d,Re1,dx,wPatch):
    
    #Assignement of coefficients on no slip walls
    if wPatch==3 or wPatch==2 or wPatch==1:
        c[1] = 0
    if wPatch==0:
        
        #Assignment of varios dx
        dx_w = dx[0][0]
        dx_e = dx[1][0]
        dx_s = dx[2][1]
        dx_n = dx[3][1]
        dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
        phi_W = (dx_e)/(dx_w*(dx_e + dx_w))
        phi_P = -(dx_e + dx_w)/(dx_e*dx_w)
        phi_E = -(phi_W + phi_P)
        
        c[0] = c[0] - c[1]*(phi_W/phi_E)
        c[4] = c[4] - c[1]*(phi_P/phi_E)
        c[1] = 0   
    return a,b,c,d

def UySouthFunction(a,b,c,d,Re1,dx,wPatch):
    #Function to assigne coefficients to points wich are adjacent to a western 
    #boundary for different boundary types
    A = np.zeros(a.shape)
    C = np.zeros(c.shape)
    D = np.zeros(d.shape)
    E = np.zeros(d.shape)
    
    #Assignment of varios dx
    dx_w = dx[0][0]
    dx_e = dx[1][0]
    dx_s = dx[2][1]
    dx_n = dx[3][1]
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    #Assignment of coefficients in the pressure boundary condition at no slip
    #walls
    A[3] = -dx_s/((dx_s + dx_n)*dx_n)
    A[4] = (dx_s + dx_n)/(dx_s*dx_n)
    A[2] = -(A[3] + A[4])
    C[3] = -2/(Re1*(dx_n + dx_s)*dx_n)
    C[4] = 2/(Re1*dx_s*dx_n)
    D[2] = 1
    E[3] = (dx_s**2)/(dx_s**2 - (dx_s + dx_n)**2)
    E[4] = -((dx_s + dx_n)**2)/(dx_s**2 - (dx_s + dx_n)**2)
    
    #Assignement of coefficients on no slip walls
    if wPatch==3 or wPatch==0: # Isothermal bottom wall or isothermal free 
                               # atmosphere
        a[4] = -(a[2]/A[2])*A[4]
        a[3] = a[3] - (a[2]/A[2])*A[3]
        c[2] = 0
        c[3] = c[3] - (a[2]/A[2])*C[3]
        c[4] = c[4] - (a[2]/A[2])*C[4]
        a[2] = 0
      
    if wPatch==2 or wPatch==1: #Adiabatic walls
        a[4] = -(a[2]/A[2])*A[4]
        a[3] = a[3] - (a[2]/A[2])*A[3]
        c[2] = 0
        c[3] = c[3] - (a[2]/A[2])*C[3]
        c[4] = c[4] - (a[2]/A[2])*C[4]
        d[3] = -(a[2]/A[2])*D[2]*E[3]
        d[4] = d[4] - (a[2]/A[2])*D[2]*E[4]
        a[2] = 0
        
    return a,b,c,d

def UyNorthFunction(a,b,c,d,Re1,dx,wPatch):
    #Function to assigne coefficients to points wich are adjacent to a western 
    #boundary for different boundary types
    A = np.zeros(a.shape)
    C = np.zeros(c.shape)
    D = np.zeros(d.shape)
    E = np.zeros(d.shape)
    
    #Assignment of varios dx
    dx_w = dx[0][0]
    dx_e = dx[1][0]
    dx_s = dx[2][1]
    dx_n = dx[3][1]
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    #Assignment of coefficients in the pressure boundary condition at no slip
    #walls
    A[2] = dx_n/((dx_s + dx_n)*dx_s)
    A[4] = -(dx_s + dx_n)/(dx_s*dx_n)
    A[3] = -(A[2] + A[4])
    C[2] = 2/(Re1*(dx_n + dx_s)*dx_s)
    C[4] = -2/(Re1*dx_s*dx_n)
    D[3] = 1
    #E[3] = (dx_s**2)/(dx_s**2 - (dx_s + dx_n)**2)
    #E[4] = -((dx_s + dx_n)**2)/(dx_s**2 - (dx_s + dx_n)**2)
    
    #Assignement of coefficients on no slip walls
    if wPatch==3 or wPatch==0: # Isothermal bottom wall or isothermal free 
                               # atmosphere
        a[4] = -(a[3]/A[3])*A[4]
        a[2] = a[2] - (a[3]/A[3])*A[2]
        c[3] = 0
        c[2] = c[2] - (a[3]/A[3])*C[2]
        c[4] = c[4] - (a[3]/A[3])*C[4]
        a[3] = 0
      
    #if wPatch==2 or wPatch==1: #Adiabatic walls
        #a[4] = -(a[2]/A[2])*A[4]
        #a[3] = a[3] - (a[2]/A[2])*A[3]
        #c[2] = 0
        #c[3] = c[3] - (a[2]/A[2])*C[3]
        #c[4] = c[4] - (a[2]/A[2])*C[4]
        #d[3] = -(a[2]/A[2])*D[2]*E[3]
        #d[4] = d[4] - (a[2]/A[2])*D[2]*E[4]
        #a[2] = 0
        
    return a,b,c,d

def Uy_dash_coefficients(U,dx,Re1,neighbours,index,wesn,bPatches):
    directions = 4
    a = np.zeros((directions+1,1))#Coefficients:p in wesnp order
    b = np.zeros((directions+1,1))#Coefficients:Ux
    c = np.zeros((directions+1,1))#Coefficients:Uy
    d = np.zeros((directions+1,1))#Coefficients:theta
    
    deltaX = dx[0][0]+dx[1][0]
    deltaY = dx[2][1]+dx[3][1]
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    #Assigning coefficients for interior points without considering boundary
    #points
    a[2] = 1/deltaY
    a[3] = -a[2]
    
    b[4] = - (U[neighbours[1]][1]-U[neighbours[0]][1])/deltaX
    
    c[0] = 2/(Re1*dx_w*deltaX) + U[index][0]/deltaX
    c[1] = 2/(Re1*dx_e*deltaX) - U[index][0]/deltaX
    c[2] = 2/(Re1*dx_s*deltaY) + U[index][1]/deltaY
    c[3] = 2/(Re1*dx_n*deltaY) - U[index][1]/deltaY
    c[4] = -2/(Re1*dx_e*dx_w) - 2/(Re1*dx_s*dx_n) - (U[neighbours[3]][1]-U[neighbours[2]][1])/deltaY
    
    d[4] = 1
    
    #Assigning coefficients considering boundary points
    if wesn[0]=='w':
        a,b,c,d = UyWestFunction(a,b,c,d,Re1,dx,bPatches[0])
    if wesn[1]=='e':
        a,b,c,d = UyEastFunction(a,b,c,d,Re1,dx,bPatches[1])
    if wesn[2]=='s':
        a,b,c,d = UySouthFunction(a,b,c,d,Re1,dx,bPatches[2])
    #if wesn[3]=='n':
        #a,b,c,d = UyNorthFunction(a,b,c,d,Re1,dx,bPatches[3])

    return a,b,c,d

def Uy_dash_equation(index,neighbours,p,U,theta,dx,wesn,bPatches,n,n_var,Re1,Re2):
    
    row = np.zeros((n_var*n, 1))

    a,b,c,d = Uy_dash_coefficients(U,dx,Re1,neighbours,index,wesn,bPatches)
    
    k=0
    for i in np.append(neighbours, index):#Among all the neighbours,
        if i<n: #if the neighbour is a inerior point
            row[i*n_var + 0] = a[k]
            row[i*n_var + 1] = b[k]
            row[i*n_var + 2] = c[k]
            row[i*n_var + 3] = d[k]
        k = k + 1
    
    return row

def thetaWestFunction(a,b,c,d,Re2,dx,wPatch):
    #Function to assigne coefficients to points wich are adjacent to a western 
    #boundary for different boundary types
    A = np.zeros(a.shape)
    C = np.zeros(c.shape)
    D = np.zeros(d.shape)
    E = np.zeros(d.shape)
    
    #Assignment of varios dx
    dx_w = dx[0][0]
    dx_e = dx[1][0]
    dx_s = dx[2][1]
    dx_n = dx[3][1]
    
    deltaX = dx[0][0]+dx[1][0]
    deltaY = dx[2][1]+dx[3][1]
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    #Assignment of coefficients in the pressure boundary condition at no slip
    #walls
    phi_P = (deltaX)/(dx_e*dx_w)
    phi_E = -(dx_w)/(dx_e*deltaX)
    
    E[4] = (phi_P)/(phi_P + phi_E)
    E[1] = (phi_E)/(phi_P + phi_E)
    
    #Assignement of coefficients on no slip walls
    if wPatch==3: # Isothermal bottom wall
        d[0] = 0
    
    if wPatch==2 or wPatch==1: #Adiabatic walls
        d[4] = d[0]*E[4] + d[4]
        d[1] = d[0]*E[1] + d[1]
        d[0] = 0
        
    if wPatch==0: #Free slip walls
        d[0] = 0
        
    return a,b,c,d

def thetaEastFunction(a,b,c,d,Re2,dx,wPatch):
    #Function to assigne coefficients to points wich are adjacent to a western 
    #boundary for different boundary types
    A = np.zeros(a.shape)
    C = np.zeros(c.shape)
    D = np.zeros(d.shape)
    E = np.zeros(d.shape)
    
    #Assignment of varios dx
    dx_w = dx[0][0]
    dx_e = dx[1][0]
    dx_s = dx[2][1]
    dx_n = dx[3][1]
    
    deltaX = dx[0][0]+dx[1][0]
    deltaY = dx[2][1]+dx[3][1]
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    #Assignment of coefficients in the pressure boundary condition at no slip
    #walls
    phi_P = -(deltaX)/(dx_e*dx_w)
    phi_W = (dx_e)/(dx_w*deltaX)
    
    E[4] = (phi_P)/(phi_P + phi_W)
    E[0] = (phi_W)/(phi_P + phi_W)
    
    #Assignement of coefficients on no slip walls
    if wPatch==3: # Isothermal bottom wall
        d[1] = 0
    
    if wPatch==2 or wPatch==1: #Adiabatic walls
        d[4] = d[1]*E[4] + d[4]
        d[0] = d[1]*E[0] + d[0]
        d[1] = 0
        
    if wPatch==0: #Free slip walls
        d[1] = 0
        
    return a,b,c,d

def thetaSouthFunction(a,b,c,d,Re2,dx,wPatch):
    #Function to assigne coefficients to points wich are adjacent to a western 
    #boundary for different boundary types
    A = np.zeros(a.shape)
    C = np.zeros(c.shape)
    D = np.zeros(d.shape)
    E = np.zeros(d.shape)
    
    #Assignment of varios dx
    dx_w = dx[0][0]
    dx_e = dx[1][0]
    dx_s = dx[2][1]
    dx_n = dx[3][1]
    
    deltaX = dx[0][0]+dx[1][0]
    deltaY = dx[2][1]+dx[3][1]
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    #Assignment of coefficients in the pressure boundary condition at no slip
    #walls
    phi_P = (deltaY)/(dx_s*dx_n)
    phi_N = -(dx_s)/(dx_n*deltaY)
    
    E[4] = (phi_P)/(phi_P + phi_N)
    E[3] = (phi_N)/(phi_P + phi_N)
    
    #Assignement of coefficients on no slip walls
    if wPatch==3: # Isothermal bottom wall
        d[2] = 0
    
    if wPatch==2 or wPatch==1: #Adiabatic walls
        d[4] = d[2]*E[4] + d[4]
        d[3] = d[2]*E[3] + d[3]
        d[2] = 0
        
    if wPatch==0: #Free slip walls
        d[2] = 0
        
    return a,b,c,d

def thetaNorthFunction(a,b,c,d,Re2,dx,wPatch):
    #Function to assigne coefficients to points wich are adjacent to a western 
    #boundary for different boundary types
    A = np.zeros(a.shape)
    C = np.zeros(c.shape)
    D = np.zeros(d.shape)
    E = np.zeros(d.shape)
    
    #Assignment of varios dx
    dx_w = dx[0][0]
    dx_e = dx[1][0]
    dx_s = dx[2][1]
    dx_n = dx[3][1]
    
    deltaX = dx[0][0]+dx[1][0]
    deltaY = dx[2][1]+dx[3][1]
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    #Assignment of coefficients in the pressure boundary condition at no slip
    #walls
    phi_P = -(deltaY)/(dx_s*dx_n)
    phi_S = (dx_n)/(dx_s*deltaY)
    
    E[4] = (phi_P)/(phi_P + phi_S)
    E[2] = (phi_S)/(phi_P + phi_S)
    
    #Assignement of coefficients on no slip walls
    if wPatch==3: # Isothermal bottom wall
        d[3] = 0
    
    if wPatch==2 or wPatch==1: #Adiabatic walls
        d[4] = d[3]*E[4] + d[4]
        d[2] = d[3]*E[2] + d[2]
        d[3] = 0
        
    if wPatch==0: #Free slip walls
        d[3] = 0
        
    return a,b,c,d

def theta_dash_coefficients(U,theta,dx,Re2,neighbours,index,wesn,bPatches):
    directions = 4
    a = np.zeros((directions+1,1))#Coefficients:p in wesnp order
    b = np.zeros((directions+1,1))#Coefficients:Ux
    c = np.zeros((directions+1,1))#Coefficients:Uy
    d = np.zeros((directions+1,1))#Coefficients:theta
    
    #Assignment of various dx
    dx_w = dx[0][0]
    dx_e = dx[1][0]
    dx_s = dx[2][1]
    dx_n = dx[3][1]
    
    deltaX = dx[0][0]+dx[1][0]
    deltaY = dx[2][1]+dx[3][1]
    dx_w, dx_e, dx_s, dx_n, deltaX, deltaY = dxCall(dx)
    
    theta_W = theta[neighbours[0]]
    theta_E = theta[neighbours[1]]
    theta_S = theta[neighbours[2]]
    theta_N = theta[neighbours[3]]
    
    U_x_P = U[index][0]
    U_y_P = U[index][1]
    
    phi = np.zeros((2,2))
    psi = np.zeros((2,3))
    
    #Assigning coefficients for interior points without considering boundary
    #points
    
    phi[0][0] = -1/deltaX
    phi[0][1] = 1/deltaX
    phi[1][0] = -1/deltaY
    phi[1][1] = 1/deltaY
    
    phi_x_W = phi[0][0]
    phi_x_E = phi[0][1] 
    phi_y_S = phi[1][0] 
    phi_y_N = phi[1][1] 
    
    psi[0][0] = 2/(dx_w*deltaX)
    psi[0][1] = -2/(dx_w*dx_e)
    psi[0][2] = 2/(dx_e*deltaX)
    psi[1][0] = 2/(dx_s*deltaY)
    psi[1][1] = -2/(dx_s*dx_n)
    psi[1][2] = 2/(dx_n*deltaY)
    
    psi_x_W = psi[0][0]
    psi_x_P = psi[0][1]
    psi_x_E = psi[0][2]
    psi_y_S = psi[1][0] 
    psi_y_P = psi[1][1] 
    psi_y_N = psi[1][2]  
    
    b[4] = -(phi_x_W*theta_W + phi_x_E*theta_E)
    
    c[4] = -(phi_y_S*theta_S + phi_y_N*theta_N)
    
    d[0] = (1/Re2)*psi_x_W - phi_x_W*U_x_P
    d[1] = (1/Re2)*psi_x_E - phi_x_E*U_x_P
    d[2] = (1/Re2)*psi_y_S - phi_y_S*U_y_P
    d[3] = (1/Re2)*psi_y_N - phi_y_N*U_y_P
    d[4] = (1/Re2)*(psi_x_P + psi_y_P)
    
    #Assigning coefficients considering boundary points
    if wesn[0]=='w':
        a,b,c,d = thetaWestFunction(a,b,c,d,Re2,dx,bPatches[0])
    if wesn[1]=='e':
        a,b,c,d = thetaEastFunction(a,b,c,d,Re2,dx,bPatches[1])
    if wesn[2]=='s':
        a,b,c,d = thetaSouthFunction(a,b,c,d,Re2,dx,bPatches[2])
    if wesn[3]=='n':
        a,b,c,d = thetaNorthFunction(a,b,c,d,Re2,dx,bPatches[3])

    return a,b,c,d

def theta_dash_equation(index,neighbours,p,U,theta,dx,wesn,bPatches,n,n_var,Re1,Re2):
    
    row = np.zeros((n_var*n, 1))

    a,b,c,d = theta_dash_coefficients(U,theta,dx,Re2,neighbours,index,wesn,bPatches)
    
    k=0
    for i in np.append(neighbours, index):#Among all the neighbours,
        if i<n: #if the neighbour is a inerior point
            row[i*n_var + 0] = a[k]
            row[i*n_var + 1] = b[k]
            row[i*n_var + 2] = c[k]
            row[i*n_var + 3] = d[k]
        k = k + 1
    
    return row