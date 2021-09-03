

**TWISTER** is a python package that helps you find and construct commensurate moiré superlattices on introducing a twist between 2D materials. Please see below for instructions regarding specific lattices. 

**TWISTER** can also be used to study the structural reconstruction of moiré superlattices with classical forcefields through a LAMMPS interface.
See FFRelaxation/ for examples

If you use this package please cite the following papers for which this code was developed:

"Twister: Construction and structural relaxation of commensurate moiré superlattices", [arxiv](https://arxiv.org/abs/2102.07884)
"Ultraflatbands and Shear Solitons in Moire Patterns of Twisted Bilayer Transition Metal Dichalcogenides", [Phys. Rev. Lett. 121, 266401 (2018)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.266401) 

We also distribute registry-dependent Kolmogorov-Crespi interlayer potentials and LAMMPS input files in the KC_ilp/ folder used in the following paper:
"Kolmogorov-Crespi Potential For Multilayer Transition Metal Dichalcogenides: Capturing Structural Transformations In Moiré Superlattices", [J. Phys. Chem. C 2019, 123, 15, 9770–9778](https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.8b10392)


For assistance with running and installation please visit [our Google Group](https://groups.google.com/g/twister-help)


## DEPENDENCIES:

python with Numpy, Scipy and Matplotlib


## HEXAGONAL LATTICES:

(With the same lattice parameter, eg. twisted bilayer graphene, 
twisted bilayer MoS2)

An analytical expression for the commensurate twist-angle 
is used. See examples/homobilayer_hex/


## ORTHORHOMBIC LATTICES and HETEROSTRUCTURES

It is often difficult to obtain exact coincidence in these systems
since lattice parameters of the constituent layers can be different. 
The commensurate twist-angles can be found using coincidence site lattice 
theory by allowing for a small strain in the constituent layers. 
See examples/heterobilayer, examples/twistedbP

