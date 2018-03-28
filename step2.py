#step2.py
# generate box and lammps input file
from class_data import data
from class_packmolgen import packmolgen
import subprocess

# prepare packmol input file
# NOTE: set box dimensions here.
box = [400., 400., 400.]
# NOTE: components and numbers
all_molecules = [   {'chain': 'chain1.xyz','bond': 'bond.x', 'angle': 'angle.x', 'Natom' : 3 , 'dihedral': 'dihedral.x', 'improper': 'improper.x' },  \
                    { 'chain': 'OH.xyz', 'Natom':48} ,\
                    {'chain':'W.xyz', 'Natom': 10 }   ]


# packmolAEM.inp is the packmol input
pack = packmolgen( boxsize = box, molecules = all_molecules, filetype = 'xyz' )

pack.packmolinpgen('MixtureAEM.xyz')

# run packmol
# NOTE: modify the packmol directory here. packmol.inp is the default
subprocess.call('/opt/packmol/packmol<packmol.inp', shell = True)

# read totalbeadlist
# necessary to convert xyz atom name into number
from necessary_beads import *
totalbeadlist = [b1, b2, b3, b4, b5, b6, b7]

# initialize d1
# (output filename, xyz of mixture, xyz of single chain, bond per chain, angle per chain, packmol input)
d1 = data('SEBS1TMA_DPD.data')

d1.load_mixture('MixtureAEM.xyz')

d1.load_all_molecules(*all_molecules)

# update box size of the lammps input file

d1.xlo = 0.
d1.xhi = box[0]
d1.ylo = 0.
d1.yhi = box[1]
d1.zlo = 0.
d1.zhi = box[2]

# NOTE: set the types manually
d1.tatom = 7
d1.tbond = 6
d1.tangle = 6
d1.tdihedral = 4
d1.timproper = 1

# NOTE
# update the total numbers of bond, angle, dihedral and improper
# must do
d1.update_totalnumber()

# write sections of lammps data file
d1.writehead()

d1.writeatom(totalbeadlist)

d1.writebond()
d1.writeangle()
d1.writedihedral()
d1.writeimproper()
