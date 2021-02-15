import sys, string
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines

def gen_pos(basis_file,alat,a1,a2,a3,sc):
  # Generates the super cell.
#  sc_pos, sc_size, n_basis = gen_supercell(basis_file,sc)
  sc_pos, sc_size, n_basis = gen_supercell_cen(basis_file,sc)
  nat= len(sc_pos)
  alat_sc =alat*sc_size
  A_ang = np.array([a1*alat_sc[0],a2*alat_sc[1],a3*alat_sc[2]])
  pos_ang, labels = crys2ang(A_ang, sc_pos,nat)
  return pos_ang,labels,nat,alat_sc,A_ang,n_basis


def locate_points_acc(a1_n,a2_n,layer1,layer2,labels1, labels2,nat,alat, file_name,celldm1):
  add = a1_n + a2_n
  delta1 = alat*a1_n/np.linalg.norm(a1_n)
  delta2 = alat*a2_n/np.linalg.norm(a2_n)
  a1_npd = a1_n + delta1 - delta2
  a1_nmd = a1_n - delta1 + delta2
  a2_npd = a2_n + delta2 - delta1
  a2_nmd = a2_n - delta2 + delta1
  add_pd = a1_n + delta1 + a2_n + delta2
  add_md = a1_n - delta1 + a2_n - delta2
  m1m2 = -1*delta1 -1*delta2
  p1p2 = delta1 + delta2
  # Create path for an outer box: V
  V = np.array([[m1m2[0],m1m2[1]],[a1_npd[0],a1_npd[1]],[add_pd[0],add_pd[1]],[a2_npd[0],a2_npd[1]]])
  bbPath = mpath.Path(V)
  # Create path for an inner box: V_in
  V_in = np.array([[p1p2[0],p1p2[1]],[a1_nmd[0],a1_nmd[1]],[add_md[0],add_md[1]],[a2_nmd[0],a2_nmd[1]]])
  bbPath_in = mpath.Path(V_in)
  encl_points = []
  pos_in_layer1 = []
  pos_in_layer2 = []
  pos_out_layer1 = []
  pos_out_layer2 = []
  fp = open(file_name,'w')
  fp.write("Superlattice vectors (celldm1_l units):\n")
  fp.write("%3.11f %3.11f %3.11f\n"%(a1_n[0]/celldm1,a1_n[1]/celldm1,a1_n[2]/celldm1))
  fp.write("%3.11f %3.11f %3.11f\n"%(a2_n[0]/celldm1,a2_n[1]/celldm1, a2_n[2]/celldm1))
  fp.write("Layer 1 points (Angstrom):\n")
  grid = bbPath.contains_points(layer1[:,0:2])
  grid_in = bbPath_in.contains_points(layer1[:,0:2])
  num = 0

  # Change in V_1.1:
  # Check overlap only for points outside inner and inside the outer box.

  for i in range(nat):
#    if bbPath.contains_point((layer1[i][0],layer1[i][1])):
    if grid[i]:
      pos = np.array([layer1[i][0],layer1[i][1],layer1[i][2]])
      Flag = False
      if grid_in[i]:
        Flag = False
      else:
        for j in range(len(pos_in_layer1)):
          pos_old = np.array([pos_in_layer1[j][1],pos_in_layer1[j][2],pos_in_layer1[j][3]])
          fold1 = pos - a1_n
          fold2 = pos - a2_n
          fold3 = pos - a2_n - a1_n
          fold4 = pos + a1_n
          fold5 = pos + a2_n
          fold6 = pos + a2_n + a1_n
          fold7 = pos + a2_n - a1_n
          fold8 = pos - a2_n + a1_n
          disp1 = np.linalg.norm(pos_old - fold1)
          disp2 = np.linalg.norm(pos_old - fold2)
          disp3 = np.linalg.norm(pos_old - fold3)
          disp4 = np.linalg.norm(pos_old - fold4)
          disp5 = np.linalg.norm(pos_old - fold5)
          disp6 = np.linalg.norm(pos_old - fold6)
          disp7 = np.linalg.norm(pos_old - fold7)
          disp8 = np.linalg.norm(pos_old - fold8)
          if disp1 < 0.5 or disp2 < 0.5 or disp3 < 0.5 or disp4 < 0.5 or disp5 < 0.5 or disp6 < 0.5 or \
             disp7 < 0.5 or disp8 < 0.5:
            Flag = True
            break
      if Flag == False:  
        fp.write("%s %3.11f %3.11f %3.11f\n"%(labels1[i],layer1[i][0],layer1[i][1],layer1[i][2]))
        pos_in_layer1.append([labels1[i],layer1[i][0],layer1[i][1],layer1[i][2]])
        num = num + 1

  n_layer1 = num
  fp.write("Number of points in layer 1: %d\n"%(n_layer1))
  fp.write("Layer 2 points (Angstrom):\n")
  grid = bbPath.contains_points(layer2[:,0:2])
  grid_in = bbPath_in.contains_points(layer2[:,0:2])
  for i in range(nat):
#    if bbPath.contains_point((layer2[i][0],layer2[i][1])):
    if grid[i]:
      pos = np.array([layer2[i][0],layer2[i][1],layer2[i][2]])
      Flag = False
      if grid_in[i]:
       Flag = False
      else:
        for j in range(len(pos_in_layer2)):
          pos_old = np.array([pos_in_layer2[j][1],pos_in_layer2[j][2],pos_in_layer2[j][3]])
          fold1 = pos - a1_n
          fold2 = pos - a2_n
          fold3 = pos - a2_n - a1_n
          fold4 = pos + a1_n
          fold5 = pos + a2_n
          fold6 = pos + a2_n + a1_n
          fold7 = pos + a2_n - a1_n
          fold8 = pos - a2_n + a1_n
          disp1 = np.linalg.norm(pos_old - fold1)
          disp2 = np.linalg.norm(pos_old - fold2)
          disp3 = np.linalg.norm(pos_old - fold3)
          disp4 = np.linalg.norm(pos_old - fold4)
          disp5 = np.linalg.norm(pos_old - fold5)
          disp6 = np.linalg.norm(pos_old - fold6)
          disp7 = np.linalg.norm(pos_old - fold7)
          disp8 = np.linalg.norm(pos_old - fold8)
          if disp1 < 0.5 or disp2 < 0.5 or disp3 < 0.5 or disp4 < 0.5 or disp5 < 0.5 or disp6 < 0.5 or \
             disp7 < 0.5 or disp8 < 0.5:
            Flag = True
            break
      if Flag == False:
        fp.write("%s %3.11f %3.11f %3.11f\n"%(labels2[i],layer2[i][0],layer2[i][1],layer2[i][2]))
        pos_in_layer2.append([labels2[i],layer2[i][0],layer2[i][1],layer2[i][2]])
        num = num + 1
  n_layer2 = num - n_layer1
  fp.write("Number of points in layer 2: %d\n"%(n_layer2))
  fp.write("\nTotal number of atoms: %d"%(num))
  print("Total number of atoms in superlattice: %d"%(num))
  print("Atom positions written to file: %s"%(file_name))
  return n_layer1, n_layer2,V, V_in

def crys2ang(A_bohr, pos,nat):
  pos_bohr = []
  labels = []
  for i in range(nat):
    pos_bohr.append(pos[i][1]*A_bohr[0] + pos[i][2]*A_bohr[1] + pos[i][3]*A_bohr[2])
    labels.append(pos[i][0])
  return np.array(pos_bohr), np.array(labels)

def read_input(input_file):
  fp = open(input_file,'r')
  lines = fp.readlines()
  fp.close()
  l2_from_file = False
  plot_lattice = True
  l2_file = ""
  for i in range(len(lines)):
    if "celldm1_l" in lines[i]:
      w = lines[i+1].split()
      alat_l = [eval(w[0]),eval(w[1]),eval(w[2])] 
    if "celldm1_u" in lines[i]:
      w = lines[i+1].split()
      alat_u = [eval(w[0]),eval(w[1]),eval(w[2])] 
    if "angle:" in lines[i]:
      w = lines[i+1].split()
      angle = eval(w[0])
    if "translate_z:" in lines[i]:
      w = lines[i+1].split()
      t_z = eval(w[0])
    if "a1_l:" in lines[i]:
      w = lines[i+1].split()
      a1_l = [eval(w[0]),eval(w[1]),eval(w[2])]
    if "a1_u:" in lines[i]:
      w = lines[i+1].split()
      a1_u = [eval(w[0]),eval(w[1]),eval(w[2])]
    if "a2_l:" in lines[i]:
      w = lines[i+1].split()
      a2_l = [eval(w[0]),eval(w[1]),eval(w[2])]
    if "a2_u:" in lines[i]:
      w = lines[i+1].split()
      a2_u = [eval(w[0]),eval(w[1]),eval(w[2])]
    if "a3_l:" in lines[i]:
      w = lines[i+1].split()
      a3_l = [eval(w[0]),eval(w[1]),eval(w[2])]
    if "a3_u:" in lines[i]:
      w = lines[i+1].split()
      a3_u = [eval(w[0]),eval(w[1]),eval(w[2])]
    if "Superlattice1:" in lines[i]:
      w = lines[i+1].split()
      SL_a1 = [eval(w[0]),eval(w[1])]
    if "Superlattice2:" in lines[i]:
      w = lines[i+1].split()
      SL_a2 = [eval(w[0]),eval(w[1])]
    if "layer2_from_file:" in lines[i]:
      w = lines[i+1].split()
      if eval(w[0]) == True:
        l2_from_file = True
        l2_file = w[1]
    if "Plot_lattice:" in lines[i]:
      w = lines[i+1].split()
      if eval(w[0]) == False:
        plot_lattice = False
  return np.array(a1_l),np.array(a2_l),np.array(a3_l),np.array(a1_u),np.array(a2_u),np.array(a3_u),angle,np.array(alat_l),np.array(alat_u),t_z,l2_from_file,l2_file,SL_a1,SL_a2,plot_lattice
  
def gen_supercell_cen(basis_file,sc):
  # This function generates a super cell from the basis atoms in crystal 
  # coordinates.
  # basis file contains the basis atoms
  # sc[0], sc[1], sc[2] integers giving the size of super cell.
  print("Generating a %d x %d x %d supercell"%(sc[0],sc[1],sc[2]))
  fp= open(basis_file,'r')
  lines = fp.readlines()
  fp.close()
  Pos = []
  for line  in lines:
    w = line.split()
    Pos.append([w[0],eval(w[1]),eval(w[2]),eval(w[3])])
  nat = len(Pos)
  sc1 = sc[0]
  sc2 = sc[1]
  sc3 = sc[2]
  Pos_sc = []
  for i in range(-1*int(sc1/2), int(sc1/2)):
    for j in range(-1*int(sc2/2), int(sc2/2)):
      for k in range(1):
        for n in range(nat):
          Pos_sc.append([Pos[n][0],Pos[n][1]/sc1 + i*(1./sc1),Pos[n][2]/sc2 + j*(1./sc2),Pos[n][3]/sc3 + k*(1./sc3)])
  return Pos_sc, np.array(sc), nat


def gen_supercell(basis_file,sc):
  # This function generates a super cell from the basis atoms in crystal 
  # coordinates.
  # basis file contains the basis atoms
  # sc[0], sc[1], sc[2] integers giving the size of super cell.
  print("Generating a %d x %d x %d supercell"%(sc[0],sc[1],sc[2]))
  fp= open(basis_file,'r')
  lines = fp.readlines()
  fp.close()
  Pos = []
  for line  in lines:
    w = line.split()
    w = filter(bool,w)
    Pos.append([w[0],eval(w[1]),eval(w[2]),eval(w[3])])
  nat = len(Pos)
  sc1 = sc[0]
  sc2 = sc[1]
  sc3 = sc[2]
  Pos_sc = []
  for i in range(sc1):
    for j in range(sc2):
      for k in range(sc3):
        for n in range(nat):
          Pos_sc.append([Pos[n][0],Pos[n][1]/sc1 + i*(1./sc1),Pos[n][2]/sc2 + j*(1./sc2),Pos[n][3]/sc3 + k*(1./sc3)])
  return Pos_sc, np.array(sc), nat
  

def gen_supercell_file(basis_file, input_file):
  fp = open(input_file,'r')
  lines = fp.readlines()
  fp.close()
  sc = []
  alat =[]
  for i in range(len(lines)):
    if "alat" in lines[i]:
      w = lines[i+1].split()
      alat = [eval(w[0]),eval(w[1]),eval(w[2])]
    if "supercell" in lines[i]:
      w = lines[i+1].split()
      sc = [eval(w[0]),eval(w[1]),eval(w[2])]
  print("Generating a %d x %d x %d supercell"%(sc[0],sc[1],sc[2]))
  fp= open(basis_file,'r')
  lines = fp.readlines()
  fp.close()
  Pos = []
  for line  in lines:
    w = line.split()
    Pos.append([w[0],eval(w[1]),eval(w[2]),eval(w[3])])
  nat = len(Pos)
  sc1 = sc[0]
  sc2 = sc[1]
  sc3 = sc[2]
  Pos_sc = []
  for i in range(sc1):
    for j in range(sc2):
      for k in range(sc3):
        for n in range(nat):
          Pos_sc.append([Pos[n][0],Pos[n][1]/sc1 + i*(1./sc1),Pos[n][2]/sc2 + j*(1./sc2),Pos[n][3]/sc3 + k*(1./sc3)])
  return Pos_sc, sc
  
def make_plot(layer1,layer2,a1,a2,a1_n,a2_n):
  add_tmp = a1_n+a2_n
  max_x = max(0,a2_n[0],add_tmp[0],a1_n[0])
  min_x = min(0,a2_n[0],add_tmp[0],a1_n[0])
  max_y = max(0,a2_n[1],add_tmp[1],a1_n[1])
  min_y = min(0,a2_n[1],add_tmp[1],a1_n[1])
  fig, ax = plt.subplots()
  aX,aY = np.array([[0, a1[0]], [0,a1[1]]])
  line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'k')
  ax.add_line(line)
  aX,aY = np.array([[0, a2[0]], [0,a2[1]]])
  line = mlines.Line2D(aX, aY , lw=1.8, alpha=1.0, color = 'k')
  ax.add_line(line)
  for i in range(len(layer1)):
    xy = (layer1[i][0],layer1[i][1])
    if xy[0] < max_x + 10 and xy[0] > min_x - 10 and xy[1] < max_y + 10 and xy[1] > min_y - 10:
      circle = plt.Circle(xy,0.24,ec="none", color = 'r')
      ax.add_patch(circle)
  
  for i in range(len(layer2)):
    xy = (layer2[i][0],layer2[i][1])
    if xy[0] < max_x + 10 and xy[0] > min_x - 10 and xy[1] < max_y + 10 and xy[1] > min_y - 10:
      circle = plt.Circle(xy,0.18,ec="none", color = 'k',alpha=1.0)
      ax.add_patch(circle)
  aX,aY = np.array([[0, a1_n[0]], [0,a1_n[1]]])
  line = mlines.Line2D(aX, aY , lw=1.3, alpha=1.0, color = 'orange')
  ax.add_line(line)
  aX,aY = np.array([[0, a2_n[0]], [0,a2_n[1]]])
  line = mlines.Line2D(aX, aY , lw=1.3, alpha=1.0, color = 'orange')
  ax.add_line(line)
  add_tmp = a1_n+a2_n
  aX,aY = np.array([[a2_n[0], add_tmp[0]], [a2_n[1],add_tmp[1]]])
  line = mlines.Line2D(aX, aY , lw=1.3, alpha=1.0, color = 'orange')
  ax.add_line(line)
  aX,aY = np.array([[a1_n[0], add_tmp[0]], [a1_n[1],add_tmp[1]]])
  line = mlines.Line2D(aX, aY , lw=1.3, alpha=1.0, color = 'orange')
  ax.add_line(line)
  plt.xlabel(r"x ($\AA$)",fontsize = 20)
  plt.ylabel(r"y ($\AA$)",fontsize = 20)
  plt.xlim(min_x-2,max_x+2)
  plt.ylim(min_y-2,max_y+2)
  
  plt.show()

def translate_coord(pos,x_sh,y_sh):
  for i in range(len(pos)):
    pos[i][0] = pos[i][0] + x_sh
    pos[i][1] = pos[i][1] + y_sh
  return pos

def Rotate_one( a1, axis, theta):
  #   Returns the rotation matrix associated with counterclockwise rotation about
  #   a given axis by  theta in radians.
  axis = np.asarray(axis)
  theta = np.asarray(theta)
  axis = axis/np.sqrt(np.dot(axis, axis))
  a = np.cos(theta/2)
  b, c, d = -axis*np.sin(theta/2)
  aa, bb, cc, dd = a*a, b*b, c*c, d*d
  bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
  M = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                   [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                   [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
  return np.dot(M,a1)

def Rotate_atoms(layer,norm,angle):
  layer_r = []
  for i in range(len(layer)): 
    layer_r.append(Rotate_one( layer[i], norm, angle))
  return np.array(layer_r)
