#!/usr/bin/python
import sys,string
import numpy as np
from funcs import *
import time

start_time = time.time()

# Read input file from twist.inp
a1_l1,a2_l1,a3_l1,a1_l2,a2_l2,a3_l2,angle,alat_l1,alat_l2,t_z,l2_from_file,l2_file,SL_a1,\
SL_a2, plot_lattice = read_input("twist.inp")

print("Inputs read from twist.inp:")
print("Unit cell lattice vectors of lower layer (normalized):")
print(a1_l1, a2_l1, a3_l1)
print("Unit cell lattice vectors of upper layer (normalized):")
print(a1_l2, a2_l2, a3_l2)
print("Lattice parameters (in Angstrom), layer 1: %f %f %f"%(alat_l1[0],alat_l1[1],alat_l1[2]))
print("Lattice parameters (in Angstrom), layer 2: %f %f %f"%(alat_l2[0],alat_l2[1],alat_l2[2]))

print("Twist angle: %f radians"%(angle))
print("Translation of layer 2 in z-direction: %f Angstroms"%(t_z))
print("Layer 2 basis read from file: %s"%(l2_from_file))

if l2_from_file:
  print("Layer 2 basis read from: %s"%(l2_file))

print("Plot the superlattice: %s"%(plot_lattice))

#alat = alat*0.52918
#t_z = t_z*0.52918

# The superlattice vectors: 
a1_ang_l1 = a1_l1*alat_l1[0]
a2_ang_l1 = a2_l1*alat_l1[1]
a1_ang_l2 = a1_l2*alat_l2[0]
a2_ang_l2 = a2_l2*alat_l2[1]

a1_n = SL_a1[0]*a1_ang_l1 + SL_a1[1]*a2_ang_l1
a2_n = SL_a2[0]*a1_ang_l1 + SL_a2[1]*a2_ang_l1
print("Superlattice vectors:")
print(a1_n,a2_n)

# Estimate the supercell size required:
tmp = a1_n + a2_n
n_uc = int(np.linalg.norm(tmp)/np.amin(alat_l1[:2]))
n_uc = int(n_uc*4.0)

# create a n_uc x n_uc x 1 supercell.
sc = [n_uc,n_uc,1]

# Place axis of rotation to the center of the 
# generated super cell.
#axis = int(n_uc/2)*a1_ang_l1 + int(n_uc/2)*a2_ang_l1
axis = [0,0]
print("Placing axis at:")
print(axis)

# Generate atom positions in the super cell for layer 1
pos_ang_l1,labels1,nat_l1,alat_sc,A_ang,n_basis = gen_pos("basis_pos_crys",alat_l1,a1_l1,a2_l1,a3_l1,sc)

layer1 = pos_ang_l1

# Generate atom positions in the super cell for layer 2
if l2_from_file:
  pos_ang_l2,labels2,nat_l2,alat_sc,A_ang,n_basis2 = gen_pos(l2_file,alat_l2,a1_l2,a2_l2,a3_l2,sc)
  layer2 = pos_ang_l2
else:
  layer2 = np.copy(pos_ang_l1)
  labels2 = np.copy(labels1)

layer1_t = np.zeros((nat_l1,3))
layer2_t = np.zeros((nat_l2,3))

# Translate layer 1 in-plane so that rotation axis is at (0,0)
layer1_t[:,0] = layer1[:,0] - axis[0]
layer1_t[:,1] = layer1[:,1] - axis[1]
layer2_t[:,0] = layer2[:,0] - axis[0]
layer2_t[:,1] = layer2[:,1] - axis[1]
layer1_t[:,2] = layer1[:,2]

# Translate layer 2 in out-of-plane direction - by interlayer spacing provided
layer2_t[:,2] = layer2[:,2] + t_z

layer1_t = np.array(layer1_t)
layer2_t = np.array(layer2_t)

# Rotate layer 2
norm = [0.,0.,1.0]
#layer2_t = layer2_t
layer2_r = Rotate_atoms(layer2_t,norm,angle)

# PLot the lattice and superlattice vectors
if plot_lattice:
  make_plot(layer1_t,layer2_r,a1_l1,a2_l1,a1_n,a2_n)

# Locate the points inside a1_n, a2_n parallelopiped (superlattice area)
nl1, nl2,V, V_in = \
locate_points_acc(a1_n,a2_n,layer1_t,layer2_r,labels1, labels2,nat_l1,alat_l1, "superlattice.dat",alat_l1[0])

## Uncomment to plot the inner and outer boxes: V, V_in
#fig, ax = plt.subplots()
#aX,aY = np.array([[0, a1_n[0]], [0,a1_n[1]]])
#line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'k')
#ax.add_line(line)
#aX,aY = np.array([[0, a2_n[0]], [0,a2_n[1]]])
#line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'k')
#ax.add_line(line)
#print V[0][0], V[1][0],V[0][1],V[1][1]
#aX,aY = np.array([[V[0][0], V[1][0]], [V[0][1],V[1][1]]])
#line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'r')
#ax.add_line(line)
#aX,aY = np.array([[V[1][0], V[2][0]], [V[1][1],V[2][1]]])
#line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'r')
#ax.add_line(line)
#aX,aY = np.array([[V[2][0], V[3][0]], [V[2][1],V[3][1]]])
#line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'r')
#ax.add_line(line)
#aX,aY = np.array([[V_in[0][0], V_in[1][0]], [V_in[0][1],V_in[1][1]]])
#line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'g')
#ax.add_line(line)
#aX,aY = np.array([[V_in[1][0], V_in[2][0]], [V_in[1][1],V_in[2][1]]])
#line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'g')
#ax.add_line(line)
#aX,aY = np.array([[V_in[2][0], V_in[3][0]], [V_in[2][1],V_in[3][1]]])
#line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'g')
#ax.add_line(line)
#plt.show()

# Check the number of atoms in each layer.
print("Number of atoms in layer 1: %d"%(nl1) )
print("Number of atoms in layer 2: %d"%(nl2) )

# Area of superlattice and the unit cell.
Ar_sl = np.linalg.norm(np.cross(a1_n,a2_n))
Ar_uc_l1 = np.linalg.norm(np.cross(a1_ang_l1,a2_ang_l1))
Ar_uc_l2 = np.linalg.norm(np.cross(a1_ang_l2,a2_ang_l2))
print("Areas, uc l1: %f, uc l2: %f, sl: %f"%( Ar_uc_l1, Ar_uc_l2, Ar_sl))
print("Expected number of atoms in layer 1 from area ratios: %d"%( round((Ar_sl/Ar_uc_l1)*n_basis)))
print("Expected number of atoms in layer 2 from area ratios: %d"%( round((Ar_sl/Ar_uc_l2)*n_basis2)))

print("Time taken: %12.2f s"%(time.time() - start_time))
