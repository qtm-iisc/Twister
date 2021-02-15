import string 

def isfloat(string):
  '''
  checks if a given string is a
  float or not;
  returns None if not;
  returns the float if it is
  '''
  success = False
  try:
    number = float(string)
    success = True
  except ValueError:
    pass
  if success == True:
    return number

def latcon(twist_input):
  '''
  Input--
  twist_input: twist.inp file
  Output--
  returns the lattice constant
  in bohr
  '''
  f = open(twist_input, 'r')
  contents = f.readlines()
  f.close()
  for i in range(len(contents)):
    if "celldm1" in contents[i]:
      a = eval(contents[i+1].split()[0])
      # return latcon in bohr
      return a * 1.889725989

def read_file(data_file):
  '''
  Input--
  data_file: superlattice.dat 
  Output--
  writes a generic input
  which can be used in qe
  with some modifications
  '''
  f = open(data_file, 'r')
  lines = f.readlines()
  f.close()
  print("-----------")
  print("Generating the input for qe")
  for i in range(len(lines)):
    if "Total number of atoms:" in lines[i]:
      tmp = lines[i].split()
      for j in range(len(tmp)):
        if isfloat(tmp[j]) is not None:
          natom = int(isfloat(tmp[j])) # total atoms
    elif "Number of points in layer 1" in lines[i]:
      tmp = lines[i].split()
      natom_l1 = int(eval(tmp[len(tmp)-1])) # atoms in layer1
    elif "Number of points in layer 2" in lines[i]:
      tmp = lines[i].split()
      natom_l2 = int(eval(tmp[len(tmp)-1])) # atoms in layer2

  a = latcon('twist.inp')
  g = open('qe.in', 'w')
  g.write('&CONTROL\n\
  calculation  = "scf"\n\
  prefix       = "lego"\n\
  pseudo_dir   = "."\n\
  outdir       = "."\n\
  verbosity="high"\n\
/\n\
&SYSTEM\n\
  ibrav     = 0\n\
  celldm(1) = %f\n\
  nat       = %d\n\
  ntyp      = 1\n\
  vdw_corr = "grimme-d2"\n\
  ecutwfc   = 60\n\
  ecutrho = 600\n\
/\n\
&ELECTRONS\n\
  conv_thr    = 1.D-8\n\
  mixing_beta = 0.3D0\n\
/\n\
&IONS\n\
 ion_dynamics="bfgs"\n\
/\n\
\n\
CELL_PARAMETERS {alat}\n'% (a, natom))
  for i in range(len(lines)):
    if "Superlattice vectors" in lines[i]: 
      w1 = lines[i+1].split()
      g.write("%s\t%s\t%s\n" %(eval(w1[0]), eval(w1[1]), eval(w1[2])))
      w2 = lines[i+2].split()
      g.write("%s\t%s\t%s\n" %(eval(w2[0]), eval(w2[1]), eval(w2[2])))
      g.write("%s\t%s\t%s\n" %('0.0', '0.0', '10.0'))
      g.write("\nATOMIC_SPECIES\nxxxxxxx\nxxxxxxx\nxxxxxx\n\nATOMIC_POSITIONS {Angstrom}\n")

    elif "Layer 1" in lines[i]:
      for j in range(i+1, natom_l1+i+1):
        w = lines[j].split()
        g.write("%s  %s  %s  %s\n" %(w[0], eval(w[1]), eval(w[2]),eval(w[3])))

    elif "Layer 2" in lines[i]:
      for k in range(i+1, natom_l2+i+1):
        w = lines[k].split()
        g.write("%s  %s  %s  %s\n" %(w[0], eval(w[1]), eval(w[2]),eval(w[3])))
  g.write('\nK_POINTS {GAMMA}')
  print('generic input for QE: qe.in')
  print('------')

# run--
read_file('superlattice.dat')
