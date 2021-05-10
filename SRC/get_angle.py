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

#Maximum supercell size

range_nm_l=X[0][0]  #lower limit 
range_nm_u=X[0][1]  #upper limit 

#Lattice 1
alat_a=X[1] #lattice parameter(Angstrom)  
#unit cell vectors 
a1 = X[2][0:2]*(alat_a[0]) 
a2 = X[3][0:2]*alat_a[1]
a3 = X[4][0:2]*alat_a[2]


#Lattice 2
alat_b= X[5] #lattice parameter(Angstrom)
# unit cell vectors
b1=X[6][0:2]*alat_b[0]
b2=X[7][0:2]*alat_b[1]
b3=X[8][0:2]*alat_b[2]

#Angle search parameters (degrees)
theta_i=X[9][0] #lower limit 
theta_f=X[9][1] #upper limit
dth=X[9][2] #step

#Tolerated mismatch (Angstrom)
mismatch=X[10]

#If TRUE,unit vectors will be strained, if FALSE lattice parameters are starined
strain_tv=X[18][0]

#Percentage strain tolerance
strain_per=X[15]/100

#Layer(s) which is(are) to be strained. Inputs: 'Top','Bottom'or'Both'
strain_layer=X[17]

#Angle between supercell vectors
fix_ang=X[11]  #will it be fixed? True or False    
f_ang=X[12]*np.pi/180  #if it is fixed, angle in radians
th_er=0.01  #if it is fixed, percentage error tolerated in the angle between the supercell vectors 

#Number of basis atoms
nba_a=X[13] #lattice 1
nba_b=X[14] #lattice 2

#Area of solution can lie within d_a percentage of the minimum area from all the possible solutions
da=0.01

#If TRUE, solutions will be plotted
plot=X[16]


#Run

SL_minAr, Ar = clt(a1,a2,b1,b2,theta_i,theta_f,dth,mismatch,strain_per,strain_layer,fix_ang,f_ang,range_nm_l,range_nm_u,th_er,alat_a,alat_b,nba_a,nba_b,plot,strain_tv)

