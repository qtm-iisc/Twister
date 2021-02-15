# -*- coding: utf-8 -*-
"""
Created on Wed May 30 17:07:57 2018

@author: saism
"""
import matplotlib.pyplot as plt
import numpy as np
from get_angle_funcs import  *
import sys

X=read_input("get_ang.inp")

# Max Supercell size

range_nm_l=X[0][0]  #Lower limit for lattice 1
range_nm_u=X[0][1]  #Upper limit for lattice 1

# Lattice 1

# lattice parameter

alat_a=X[1]  

# unit cell vector

a1 = X[2][0:2]*(alat_a[0]) 
a2 = X[3][0:2]*alat_a[1]
a3 = X[4][0:2]*alat_a[2]

basis1=[[0,0]]


#Lattice 2

# lattice parameter

alat_b= X[5]

# unit cell vector

b1=X[6][0:2]*alat_b[0]
b2=X[7][0:2]*alat_b[1]
b3=X[8][0:2]*alat_b[2]

basis2=[[0,0]]

#Angle search 

theta_i=X[9][0] #lower limit 
theta_f=X[9][1] #upper limit
dth=X[9][2] #step

# percentage strain tolerance

mismatch=X[10]
strain_tv=X[18][0]
strain_per=X[15]/100
strain_layer=X[17]


#Angle between supercell vectors

fix_ang=X[11]  #True or False    

f_ang=X[12]*np.pi/180  #angle

th_er=0.01  #error

#Number of basis atoms

nba_a=X[13]
nba_b=X[14]

#percentage error in min area

da=0.01

#plot
plot=X[16]


#run

SL_minAr, Ar = clt(a1,a2,b1,b2,theta_i,theta_f,dth,mismatch,strain_per,strain_layer,fix_ang,f_ang,range_nm_l,range_nm_u,th_er,alat_a,alat_b,nba_a,nba_b,plot,strain_tv)

