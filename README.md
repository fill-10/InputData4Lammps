# DPDLammpsInput

This is a very simple tool to generate the lammps input data file for DPD simulations.

User defines everything. All constrains are defined manually.

Tested on python3.6. No py plug-ins are required. However, it requires Packmol.

HOWTO:

Create all atoms you want to use in necessary_beads.py

Use step1.py to create single molecules and their bond/angle/dihedral/impropler list files.

Use step2.py to generate the box by Packmol, and write the lammps data file. 

Existing problems:

Masses are defined as 1.0 in the lammps data file. Please update masses in lammps input script.

In xyz files, atom names are stored as real names but not the numbers as in lammps for the visualization purposes.

############

Related publications:

[ Please cite ]

[1] Luo X.; Paddison S.J, DPD simulations of anion exchange membrane: The effect of an alkyl spacer on the hydrated morphology, Solid State Ionics, 2019, 339, 115012

[2] Luo X.; Liu, H.; Bae, C.; Tuckerman, M.E.; Hickner, M.A.; Paddison S.J., Mesoscale Simulations of Quaternary Ammonium-Tethered Triblock Copolymers: Effects of the Degree of Functionalization and Styrene Content, The Journal of Physical Chemistry C, 2020, 124 (30), 16315-16323


X. Luo
