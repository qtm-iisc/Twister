#!/usr/bin/python
import matplotlib.lines as mlines
import numpy as np
import matplotlib.pyplot as plt
from funcs import *

def readfile(file_name):
  fp = open(file_name,'r')
  lines = fp.readlines()
  fp.close()
  A = []
  for i in range(len(lines)):
    w = lines[i].split()
    A.append([w[0], eval(w[1]),eval(w[2]),eval(w[3])])
  return A

alat = 6.2866584 
celldm1 = alat
A1_n = np.array([16.46207796116296,0.,0.])*alat
A2_n = np.array([0.,28.513154847096512,0.])*alat

a1 = np.array([0.5, 0.8660254, 0.0])*alat/2
a2 = np.array([-0.5, 0.8660254, 0.0])*alat/2
a1_n = 9*a1 + 10*a2
a2_n =-10*a1 + 19*a2
pos_l = readfile('pos_l')
pos_u = readfile('pos_u')

labels1 = []
layer1 = []
labels2 = []
layer2 = []
for p in pos_l:
  if 'Mo' in p[0]:
    labels1.append('Mo')
  elif 'S' in p[0]:
    labels1.append('S')
  layer1.append([p[1],p[2],p[3]])
 
for p in pos_u:
  if 'Mo' in p[0]:
    labels2.append('Mo')
  elif 'S' in p[0]:
    labels2.append('S')
  layer2.append([p[1],p[2],p[3]])
nat = len(layer1)
layer1 = np.array(layer1)
layer2 = np.array(layer2)

fig, ax = plt.subplots()
plt.scatter(layer1[:,0],layer1[:,1],s = 5, color =  'darkgreen')
plt.scatter(layer2[:,0],layer2[:,1],s = 5,color = 'navy')
aX,aY = np.array([[0, A1_n[0]], [0,A1_n[1]]])
line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'k')
ax.add_line(line)
aX,aY = np.array([[0, A2_n[0]], [0,A2_n[1]]])
line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'k')
ax.add_line(line)
aX,aY = np.array([[0, a1_n[0]], [0,a1_n[1]]])
line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'r')
ax.add_line(line)
aX,aY = np.array([[0, a2_n[0]], [0,a2_n[1]]])
line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'r')
ax.add_line(line)
plt.show()

locate_points_acc(A1_n,A2_n,layer1,layer2,labels1, labels2,nat,alat, 'enclosed_points.dat',celldm1)
