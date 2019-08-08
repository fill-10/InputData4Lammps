# DPDLammpsInput

This is a very simple tool to generate the lammps input data file for DPD simulations.

User defines everything. All constrains are defined manually.

It is only dependent on the standard python2 environment. Tested on python2.7.

No py plug-ins are required.

However, it requires Packmol.

HOWTO:

Create all atoms you want to use in necessary_beads.py

Use step1.py to create single molecules and their bond/angle/dihedral/impropler list files.

Use step2.py to generate the box by Packmol, and write the lammps data file. 

Existing problems:

Masses are defined as 1.0 in the lammps data file. Please update masses in lammps input script.

In xyz files, atom names are stored as real names but not the numbers as in lammps.

This is for the visualization purposes.

############

Please cite: 
Luo X.; Paddison S.J, DPD simulations of anion exchange membrane: The effect of an alkyl spacer on the hydrated morphology, Solid State Ionics 339, 115012

X. Luo
