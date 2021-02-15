#!/usr/bin/python

import numpy as np

fp = open("twist.inp", 'r')
lines= fp.readlines()
fp.close()

for i in range(len(lines)):
  if "celldm1_l" in lines[i]:
    w = lines[i+1].split()
    alat = eval(w[0]) 

fp = open("superlattice.dat", 'r')
lines= fp.readlines()
fp.close()


w= lines[1].split()
A1 = np.array([eval(w[0]), eval(w[1]), eval(w[2])])*alat
w= lines[2].split()
A2 = np.array([eval(w[0]), eval(w[1]), eval(w[2])])*alat
A3 = np.array([0.0, 0.0, 25.0])


for i in range(len(lines)):
  if "Number of points in layer 1" in lines[i]:
    w = lines[i].split()
    n_l1 = eval(w[6])
  if "Number of points in layer 2" in lines[i]:
    w = lines[i].split()
    n_l2 = eval(w[6])

spcs = []
n_spc = []
pos = []
for i in range(4, 4 +n_l1):
  w = lines[i].split()
  if w[0] in spcs:
    pos.append([w[0], eval(w[1]), eval(w[2]), eval(w[3])])
    n_spc[spcs.index(w[0])] += 1 
  else:
    spcs.append(w[0])
    pos.append([w[0], eval(w[1]), eval(w[2]), eval(w[3])])
    n_spc.append(1)

for i in range(6 + n_l1, 6 +n_l1 + n_l2):
  w = lines[i].split()
  if w[0] in spcs:
    pos.append([w[0], eval(w[1]), eval(w[2]), eval(w[3])])
    n_spc[spcs.index(w[0])] += 1 
  else:
    spcs.append(w[0])
    pos.append([w[0], eval(w[1]), eval(w[2]), eval(w[3])])
    n_spc.append(1)



print("Species,  number of atoms: ", spcs, len(pos))

fp = open("POSCAR.vasp",'w')
fp.write("Label\n")
fp.write("1.0\n")
fp.write("%3.12f %3.12f %3.12f\n"%(A1[0], A1[1],A1[2]))
fp.write("%3.12f %3.12f %3.12f\n"%(A2[0], A2[1],A2[2]))
fp.write("%3.12f %3.12f %3.12f\n"%(A3[0], A3[1],A3[2]))
for i in range(len(spcs)):
  fp.write("%s "%(spcs[i]))
fp.write("\n")
for i in range(len(spcs)):
  fp.write("%d "%(n_spc[i]))
fp.write("\n")

fp.write("Cartesian\n")
for sp in spcs:
  for i in range(n_l1 + n_l2):
    if pos[i][0] == sp:
      fp.write("%3.12f %3.12f %3.12f\n"%(pos[i][1], pos[i][2], pos[i][3]))

fp.close()
