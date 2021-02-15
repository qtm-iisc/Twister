# lammps to perform relaxation
from lammps import lammps
lmp = lammps()

# run lammps
lmp.file("lammps.in")
