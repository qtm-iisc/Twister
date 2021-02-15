import numpy as np

f = open("lammps.dat_min", "r")
lines = f.readlines()
f.close()

atom_id = []
atom_type = []
x = []
y = []
z = []

g = open("relaxed_data", "w")
for i in range(len(lines)):
  if "atoms" in lines[i]:
    natom = int(eval(lines[i].split()[0]))
    g.write("%s\n%d\n"%("natom :", natom))
  elif "atom types" in lines[i]:
    at_type = int(eval(lines[i].split()[0]))
  elif "xlo xhi" in lines[i]:
    xlo = eval(lines[i].split()[0])
    xhi = eval(lines[i].split()[1])
    g.write("\nMoire lattice vectors (in Angstrom)\n")
    g.write("%s\n%.10f %.10f %.10f\n"%("a1 :", (xhi-xlo), 0.0, 0.0))
  elif "ylo yhi" in lines[i]:
    ylo = eval(lines[i].split()[0])
    yhi = eval(lines[i].split()[1])
  elif "xy xz yz" in lines[i]:
    xy = eval(lines[i].split()[0])
    xz = eval(lines[i].split()[1])
    yz = eval(lines[i].split()[2])
    g.write("%s\n%.10f %.10f %.10f\n"%("a2 :", xy, yhi-ylo, 0.0))
  elif "zlo zhi" in lines[i]:
    zlo = eval(lines[i].split()[0])
    zhi = eval(lines[i].split()[1])
    g.write("%s\n%.10f %.10f %.10f\n"%("a3 :", 0.0, 0.0, zhi-zlo)) # 2D
  elif "Masses" in lines[i]:
    g.write("\n%s\n"%("type    Mass"))
    for j in range(i+2, at_type+i+2):
      g.write("%d\t%.10f\n"%(int(eval(lines[j].split()[0])), eval(lines[j].split()[1])))
  elif "Atoms" in lines[i] or \
    "Atoms # atomic" in lines[i] :
    g.write("\n%s\n%s\n"%("Atomic positions (in Angstrom)", "id    type    x    y    z"))
    for j in range(i+2, natom+i+2):
      g.write("%d\t%d\t%.10f\t%.10f\t%.10f\n"%(int(eval(lines[j].split()[0])), int(eval(lines[j].split()[1])), eval(lines[j].split()[2]), eval(lines[j].split()[3]), eval(lines[j].split()[4])))
    
    
g.close()
