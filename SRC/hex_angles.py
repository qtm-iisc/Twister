#!/usr/bin/python

import sys, string
import numpy as np
alat = 2.46 # Angstroms
i_ang = []
for i in range(41):
  angle =  np.arccos((3*i**2 + 3*i + 0.5)/(3*i**2 + 3*i + 1))*180/np.pi
  print i, angle,  (i,i+1), (-1*(i+1),2*i + 1), angle*np.pi/180, np.sqrt(3*i**2 + 3*i + 1)*alat
