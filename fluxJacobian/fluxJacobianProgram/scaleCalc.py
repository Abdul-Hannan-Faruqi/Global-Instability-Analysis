# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 20:20:05 2020

@author: Ahmad Faraz
"""
import numpy as np

#Properties of Air at 310K
alpha = 2.15e-5
nu_kinematic = 1.53e-5
g = 9.81
beta = 3.37e-3
deltaT = 20 #Temperature difference
Pr = 0.7116
rho = 1.19
k = 0.0258

AspectRatio = 1
#L_cavity = AspectRatio*H_cavity

def L_scale(Ra):
    H_cavity = ((Ra*alpha*nu_kinematic)/(g*beta*deltaT))**(1/3)
    return H_cavity

def U_scale(Ra):
    H_cavity = ((Ra*alpha*nu_kinematic)/(g*beta*deltaT))**(1/3)
    return (g*beta*H_cavity*deltaT)**(1/2)

def p_scale(Ra):
    H_cavity = ((Ra*alpha*nu_kinematic)/(g*beta*deltaT))**(1/3)
    return rho*g*beta*H_cavity*deltaT
