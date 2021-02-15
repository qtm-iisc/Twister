#!/bin/bash

# Modify the following few lines 
#------------------------------------------------
# Path to twister.py, auxillary files and basis position
mos2_mose2="/home/user/codes/GetAng/updated_examples/mos2_mose2"
twister_path="/home/user/codes/GetAng/SRC/twister"
aux_files="$mos2_mose2/aux_files"
basis_file="$aux_files/basis_0"
# to perform lammps relaxation
# set lammps_relax=false if you don't want to relax
# If set false, forcefield, style, and lammps_path are not
# used; If you don't want to use python, style=""
lammps_relax=false
forcefield="/home/user/codes/GetAng/bP"
style=""
lammps_path="/home/user/codes/lammps-30Mar18/src"
#-------------------------------------------------


# rigidly twisted structures creation
cp $basis_file/basis_pos_crys* .
python3 $twister_path/twister.py twist.inp 

# generic quantum espresso input file
cp $aux_files/toqe.py ./
python3 toqe.py
rm toqe.py

# lammps.dat file creation  
cp $aux_files/lammps_triclinic.py ./
cp $aux_files/Mass_FF ./
python3 lammps_triclinic.py
rm lammps_triclinic.py

# prepare lammps input file
cp $aux_files/tolammps.py ./
python3 tolammps.py
rm tolammps.py

# lammps minimization (serial) 
if [[ "$lammps_relax" = true ]]; then
echo "will perform lammps minimization"
cp $forcefield/* .
cp $aux_files/extract_pos.py .
cp $aux_files/run_lammps.py .
if [[ "$style" = "use_python" ]]; then
echo "will be using python wrapper for lammps" 
python3 run_lammps.py
python3 extract_pos.py
else
echo "will be using standard lammps executable"
$lammps_path/lmp_serial -in lammps.in
cp $aux_files/extract_pos.py .
python3 extract_pos.py
fi
fi
echo "==Done=="

