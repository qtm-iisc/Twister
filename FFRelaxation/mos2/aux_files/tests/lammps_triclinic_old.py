import numpy as np 
import string 

def get_names(basis_file):
  '''
  Input--
  basis_file: atom_type, positions
  in the crystal coordinate should 
  be given
  Output--
  returns the strings of atom name
  ''' 
  a, b, c, d = np.loadtxt(basis_file, unpack=True, dtype='str')
  return a

def get_mass(basis_file1, basis_file2, mass_file):
  '''
  Input--
  mass_file: files containing masses
  basis_file1, basis_file2
  Output--
  Returns the masses of different atom types
  '''
  f = open(mass_file, 'r') 
  lines = f.readlines()
  f.close()
  basis1 = get_names(basis_file1).tolist()
  basis2 = get_names(basis_file2).tolist()
  M = []
  for i in range(len(lines)):
    if "Masses_l1" in lines[i]:
      for j in range(i+1, i+1+len(np.unique(basis1))):
        M.append(eval(lines[j].split()[1]))
    elif "Masses_l2" in lines[i]:
      for j in range(i+1, i+1+len(np.unique(basis2))):
        M.append(eval(lines[j].split()[1]))
  return M
      

def read_qe(input_file, data_file):
  '''
  Input--
  input_file: qe.in,
  data_file: pos_crys_mod
  Returns--
  lattice vectors, natom, lattice 
  constant
  '''
  f = open(input_file, 'r')
  lines = f.readlines()
  f.close()
  A = []
  at_type = []
  pos = []
  for i in range(len(lines)):
    if "CELL_PARAMETERS" in lines[i] or "cell_parameters" in lines[i]:
      for j in range(3):
        w = lines[i+j+1].split()
        A.append( [eval(w[0]), eval(w[1]), eval(w[2])])
    elif "nat" in lines[i]:
      natom = eval(lines[i].split()[2])
    elif "celldm(1)" in lines[i]:
      alat = eval(lines[i].split()[2])
    elif "ATOMIC_POSITIONS" in lines[i] or "atomic_positions" in lines[i]:
      for j in range(natom):
        w = lines[i+j+1].split()
        at_type.append(w[0])
        pos.append([eval(w[1]), eval(w[2]), eval(w[3])])

      # converting backto angstrom for lammps metal unit
      A = np.array(A)*alat*0.52918
      pos = np.array(pos)
      Ainv = np.linalg.inv(A)
      g = open(data_file,'w')
      for i in range(natom):
        p_c = np.dot(pos[i], Ainv)
        # to wrap
        #g.write("%s %f %f %f\n"%(at_type[i], np.mod(p_c[0], 1.0), np.mod(p_c[1], 1.0), p_c[2]))
        # no wrap
        g.write("%s %f %f %f\n"%(at_type[i], p_c[0], p_c[1], p_c[2]))
      g.close()
  return  A, natom,  alat


def rot_matrix(input_file, data_file):
  '''
  Input--
  input_file: qe.in
  data_file: pos_crys_mod
  Returns--
  '''
  A, natom, alat = read_qe(input_file, data_file)
  a1 = A[0]/np.linalg.norm(A[0])
  a2 = A[1]/np.linalg.norm(A[1])

  # theta needs to be rotated to align 
  # a1 along the x-axis (counterclock-wise)
  # They are in radian in the following
  x_ax = np.array([1.0, 0.0, 0.0])
  if a2[1] >= 0 :
    theta = np.arccos(np.dot(a1, x_ax))  + np.pi 
  elif a2[1] < 0 :
    theta = np.arccos(np.dot(a1, x_ax))
  R = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 1], [0, 0, 1]])       # R matrix
  a1_t = np.linalg.solve(R, a1)
  a2_t = np.linalg.solve(R, a2)
  R_inv = np.linalg.inv(R)

  print("\n\n--Cross-check the rotation--\n")
  print("original a1 and transformed a1 for lammps input: ", a1.tolist(), a1_t.tolist())
  print("The transformed a1 should be along x-axis (+/-): ", x_ax.tolist())
  print("original a2 and transformed a2 for lammps input :", a2.tolist(), a2_t.tolist())
  print("\n----------------------\n\n")

  return A, np.array([a1_t*np.linalg.norm(A[0]), a2_t*np.linalg.norm(A[1]), A[2]]), natom, alat


def lammps_data(input_file, data_file):
  '''
  Input--
  input_file: qe.in, 
  data_file: pos_crys_mod
  Output--
  creates a lammps.dat file
  This is *only* targeted to 
  use for 2D materials
  '''
  # axis alignements by rotating 
  v1, v2, natom, alat = rot_matrix(input_file, data_file)
  # load data from data_file
  at_type, x_c, y_c, z_c = np.loadtxt(data_file, dtype={'names': ('label', 'x_c', 'y_c', 'z_c'), 'formats': ('|S15', np.float, np.float, np.float)}, unpack = True)
  # atom position in Angstroms
  pos_A_x = [np.dot(x_c[i]*v2[0], np.array([1.0,0.0,0.0])) + np.dot(y_c[i]*v2[1], np.array([1.0,0.0,0.0])) + np.dot(z_c[i]*v2[2], np.array([1.0,0.0,0.0])) for i in range(natom)]
  pos_A_y = [np.dot(x_c[i]*v2[0], np.array([0.0,1.0,0.0])) + np.dot(y_c[i]*v2[1], np.array([0.0,1.0,0.0])) + np.dot(z_c[i]*v2[2], np.array([0.0,1.0,0.0])) for i in range(natom)]
  pos_A_z = [np.dot(x_c[i]*v2[0], np.array([0.0,0.0,1.0])) + np.dot(y_c[i]*v2[1], np.array([0.0,0.0,1.0])) + np.dot(z_c[i]*v2[2], np.array([0.0,0.0,1.0])) for i in range(natom)]
   
  # basis names taken as strings
  basis1 = get_names('basis_pos_crys').tolist()
  basis2 = get_names('basis_pos_crys_l2').tolist()
  atom_types = len(np.unique(basis1)) + len(np.unique(basis2))
  print("")
  print("Atom types are inferred from basis file")
  print("atom_types: \n%d"%(atom_types))
  print("")

  # Masses 
  M = get_mass("basis_pos_crys", "basis_pos_crys_l2", "Mass_FF")

  # simulation box dimensions
  print("Creating the lammps data file: lammps.dat")
  # The in-plane lattice vectors
  ax1 = v2[0]/np.linalg.norm(v2[0])
  ax2 = v2[1]/np.linalg.norm(v2[1])
  # angle between ax1, ax2
  tilt_ang = np.arccos(np.dot(ax1, ax2))
  print("")
  print("\nInformation about simulating box-dimensions")
  print("Lx (xlo, xhi) : ", 0.00000, np.linalg.norm(v2[0]))	
  print("tilt (xy) : ", np.cos(tilt_ang) * np.linalg.norm(v2[1]))
  print("Ly (ylo yhi) : ", 0.0000, np.sin(tilt_ang) * np.linalg.norm(v2[1]))
  # 2D materials with a large vacuum spacing
  print("Lz (zlo zhi) :", 0.0, 100.0)
  print("See lammps documentation for details.\n")
  print("")
  
  # the lammps data file generation 
  f = open("lammps.dat",'w')
  f.write("\n %d %s"%(natom, 'atoms'))
  
  f.write("\n %d %s\n"%(atom_types, 'atom types'))
  f.write("\n %f %f %s"%(0.00000, np.linalg.norm(v2[0]), 'xlo xhi'))
  f.write("\n %f %f %s"%(0.00000, np.sin(tilt_ang)*np.linalg.norm(v2[1]), 'ylo yhi'))
  f.write("\n %f %f %s"%(0.00000, 100.0, 'zlo zhi'))
  f.write("\n %f %f  %f %s\n"%(np.cos(tilt_ang)*np.linalg.norm(v2[1]), 0.0, 0.0, 'xy xz yz'))
  f.write("\n %s \n \n"%('Masses'))
  for i in range(1, atom_types+1):
    # atom_type 1 corresponds to M[0]
    # just counting initilizations
    f.write("%d %f \n"%(i, M[i-1] ))
  f.write("\n %s \n \n"%('Atoms'))


  # to distinguish top and bottom layers
  z_mean = np.average(pos_A_z)
  for i in range(natom):
    if pos_A_z[i] < z_mean:
      indx = basis1.index(at_type.tolist()[i].decode()) + 1   # top layer atom_type decoded from bytes
    if pos_A_z[i] > z_mean: 
      indx = basis2.index(at_type.tolist()[i].decode()) + 1 + len(np.unique(at_type)) # bottom layer atom type decoded
    f.write("%d %d %f %f %f\n"%(i + 1, indx, pos_A_x[i], pos_A_y[i], pos_A_z[i]))
  print("lammps.dat file is written")
  print("--------")
  f.close()	

# ---- run -------
lammps_data('qe.in', 'pos_crys_mod')
