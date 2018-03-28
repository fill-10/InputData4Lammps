# DPDLammpsInput

This is a very simple tool to generate the lammps input data file for md/dpd simulation.

User defines everything. All constrains are defined manually.

It is only dependent on the standard python3 environment. Tested on python3.6.

No py plug-ins are required.

However, it requires Packmol.

HOWTO:

Write atoms in necessary_beads.py

Use step1.py to create single molecules and their bond/angle/dihedral/impropler list files.

Use step2.py to generate the box by Packmol, and write the lammps data file. 

############
Have fun!

X. Luo
