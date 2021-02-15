# lammps input file

f = open('lammps.in','w')
print("Preparing the input file for lammps: lammps.in")
f.write("#Initialize--\n\
#general settings\n\
units           metal\n\
dimension       3\n\
box tilt        large\n\
atom_style      atomic\n\
\n\
# structure\n\
boundary        p p p\n\
read_data       lammps.dat\n\
\n")

g = open("Mass_FF", "r")
lines=g.readlines()
g.close()
for i in range(len(lines)):
  if "FFs" in lines[i]:
    indx = i
for i in range(indx+1, len(lines)):
  f.write(lines[i])

f.write("neighbor        2.0 bin\n\
neigh_modify every 1 delay 0 check yes\n\
\n\
#optimize at 0 K\n\
dump            1 all custom 100 dump.initial id type x y z\n\
thermo          1000\n\
thermo_style    custom step pe press\n\
undump          1\n\
\n\
min_style       fire\n\
minimize        0.0 1.0e-4 1000000 1000000\n\
write_data      lammps.dat_min\n\
")
f.close()
