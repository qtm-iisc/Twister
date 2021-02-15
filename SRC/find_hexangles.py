#!/usr/bin/python

import sys, string
import numpy as np
alat = 3.14 # Angstroms
i_ang = []
for m in np.arange(1,120):
  for n in np.arange(1,120):
    m = np.float(m)
    n = np.float(n)
    angle =   np.arccos((m**2 + 4*m*n + n**2)/(2*(m**2 + m*n + n**2)))*180./np.pi
    if angle > 29 and angle < 30:
      print angle,  (m,n), (m+n,-1*m), angle*np.pi/180, np.sqrt(m**2 + m*n + n**2)*alat
