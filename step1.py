#! /usr/bin/python
# step1.py
# generate xyz, bond, angle for a single chain

from class_bead import bead
from class_monomer import monomer
from class_block import block
from class_polymer import polymer
from class_data import data
#from class_bond import Bond
#from class_angle import Angle
from class_packmolgen import packmolgen
import subprocess

# NOTE 1
# create beads
# b= bead(Name, coordinates, namenumber)

b1 = bead('S', [0., 1.0, 0.],2)

b2 = bead('B', [0., 0.0, 0.],1)

b1p = bead('S', [0.,-1.0, 0.],2)

b4 = bead('M', [0, -2.0, 0.],4)

b5 = bead('A', [0.,-3.0, 0.],5)

#

#

b3 = bead('E',[0.0, 0.0, 0.0], 3)

#

b6 = bead('O', [0., 0., 0.], 6)

b7 = bead('W', [0., 0., 0.],7)


# NOTE 2: define monomer
# m1=monomer(Name, mon index, beadthread)

beadthread = [b1, b2, b1p, b4, b5]

m1 = monomer('M1',1 ,beadthread)

# update topology. 
# monomerobject.update_first(beadindex)
# also available for last, second, lastbut1
m1.update_first(2)
m1.update_second(2)
m1.update_last(2)
m1.update_lastbut1(2)

# bond numbers are consistent to Fatemeh's note
# bondlist =[ [bondtype, bead1_index, bead2_index],... ] nomer1
m1.bondlist = [[2,1,2],[2,2,3],[5,3,4],[6,4,5]]
m1linkbond = [1, None,2]
m1.bondlist.append(m1linkbond)


# anglelist = [ [angletype, bead1, beadapex, bead3] ]
m1.anglelist = [[5, 2, 3, 4],[6, 3, 4, 5]]
m1linkangles= [[1,None,2,None]]
m1.anglelist += m1linkangles
m1.dihedrallist = [[6, -1, -2, 2, 1], [7, 1,2,3,4] ]


# NOTE 3: define block
# define linking bond of m1 to m1 inside block1 (bk1)

m1linkvector = [1.0 , 0.0, 0.0]

bk1 = block('bk1', m1, 8, m1linkvector)

#print bk1.deltaV

bk1.prebond = None

# NOTE 4: define other blocks
#########################
m2 = monomer('m2', 2, [b3])
m2.bondlist = [ [4,None,1] ]
m2.anglelist =[ [4,None, 1, None] ]

m2.dihedrallist= []

m2linkvector = [1.0, 0.0, 0.0]
bk2 = block('bk2', m2, 152, m2linkvector)
bk2.prebond = [ 3, None,1]
bk2.prelinkvector = [1.0, 0.0, 0.0]
bk2.headangle_list = [[2, None, None, 1], [3, None, 1, 1] ]
bk2.headdih_list = [[1, -2, -2, -2, 1]]

#########################
bk3 = block('bk3', m1, 8, m1linkvector)
bk3.prebond = [3, None, 2]
bk3.prelinkvector = [1.0, 0.0, 0.0]
bk3.headangle_list = [ [ 3, None, None, 2], [2, None, 2, 2] ]
bk3.headdih_list = [[2, -1, -1, 2, 1]]


# NOTE 5: define polymer chain, polymer( name, blocklist )

blocklist1 = [ bk1, bk2, bk3 ]

p1 = polymer('p1', blocklist1)


# write three files

# write the xyz file of chain
p1.writexyz('chain1.xyz')


# write the bond file
# p1.writebond(filename)

p1.writebond('bond.x')

# write the angle file
# p1.writeangle(filename)

p1.writeangle('angle.x')


p1.writedihedral('dihedral.x')


# prepare packmol input file
# NOTE: One can set box dimensions here,  for example.
box = [400., 400., 400.]
# NOTE: components and numbers
components=[['chain1.xyz',3], ['OH.xyz', 48], ['W.xyz', 10]]

# packmolAEM.inp is the packmol input
pack = packmolgen('packmolAEM.inp', box, components)

pack.packmolinpgen('MixtureAEM.xyz')

# run packmol
# NOTE: modify the packmol directory here
subprocess.call('/opt/packmol/packmol<packmolAEM.inp', shell = True)

# run step2

totalbeadlist = [b1, b2, b3, b4, b5, b6, b7]
totalanglelist = []

# initialize d1
# (output filename, xyz of mixture, xyz of single chain, bond per chain, angle per chain, packmol input)
d1 = data('SEBS1TMA_random_DPD.data')

d1.load_mixture('MixtureAEM.xyz')
d1.load_chain('chain1.xyz')
d1.load_nmol('packmolAEM.inp')
d1.load_bond('bond.x')
d1.load_angle('angle.x')
d1.load_dihedral('dihedral.x')


d1.xlo = 0.
d1.xhi = box[0]
d1.ylo = 0.
d1.yhi = box[1]
d1.zlo = 0.
d1.zhi = box[2]

# One can set types of atom, bond, angle, dihedral, improper
# NOTE
d1.tatom = 7
d1.tbond = 6
d1.tangle = 6
d1.tdihedral = 4
#d1.timproper = 0

d1.writehead()

# writeatom(beadnamelist) outputs the atom cooridnates.

d1.writeatom(totalbeadlist)

d1.writebond()
d1.writeangle()
d1.writedihedral()
