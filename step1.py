# step1.py
# generate xyz, bond, angle for a single chain

#from class_bead import bead
from class_monomer import monomer
from class_block import block
from class_polymer import polymer

# NOTE 1
# create beads
# b= bead(Name, coordinates, namenumber)
# beads in necessary_beads.py
from necessary_beads import *

# NOTE 2: define monomer
# m1=monomer(Name, mon index, beadthread)

beadthread = [b1, b2, b1p, b4, b5]

m1 = monomer(beadthread)

# must update
# the first and last
# this will determine the deltaV and writexyz
m1.update_first(2)
m1.update_last(2)

# bond numbers are consistent to Fatemeh's note
# bondlist =[ [bondtype, bead1_index, bead2_index],... ] nomer1
m1.bondlist = [[2,0,1],[2,1,2],[5,2,3],[6,3,4]]
m1linkbond = [1, -4 ,1]
m1.bondlist.insert(0,m1linkbond)


# anglelist = [ [angletype, bead1, beadapex, bead3] ]
m1.anglelist = [[5, 1, 2, 3],[6, 2, 3, 4]]
m1linkangles= [1,-9,-4 ,1]
m1.anglelist.insert(0, m1linkangles)
m1.dihedrallist = [[3, -1, -2, 1, 0], [4, 0,1,2,3] ]


# NOTE 3: define block
# define linking bond of m1 to m1 inside block1 (bk1)

m1linkvector = [1.0 , 0.0, 0.0]

bk1 = block(m1, 8, m1linkvector)

#print( bk1.deltaV)

bk1.prebond = [ ] 

# NOTE 4: define other blocks
#########################
m2 = monomer([b3])
m2.bondlist = [ [4,-1 ,0] ]
m2.anglelist =[ [4,-2, -1, 0] ]

m2.dihedrallist= []

m2linkvector = [1.0, 0.0, 0.0]
bk2 = block( m2, 152, m2linkvector)
bk2.prebond = [ [ 3, -4, 0] ]
bk2.prelinkvector = [1.0, 0.0, 0.0]
bk2.headangle_list = [[2, -9, -4, 0], [3,-4 , 0, 1] ]
bk2.headdih_list = [[1, -14, -9, -2, 0]]
bk2.headimproper_list=[[1, -10, -9, -4 , 0]]

#########################
bk3 = block( m1, 8, m1linkvector)
bk3.prebond = [ [3, -1, 1] ]
bk3.prelinkvector = [1.0, 0.0, 0.0]
bk3.headangle_list = [ [ 3, -2, -1, 1], [2, -1, 1, 6] ]
bk3.headdih_list = [[2, -2, -1, 1, 0]]


# NOTE 5: define polymer chain, polymer( name, blocklist )

blocklist1 = [ bk1, bk2, bk3 ]

p1 = polymer(blocklist1)


# write the xyz file of chain
p1.writexyz('chain1.xyz')


# write the bond file
p1.writebond('bond.x')

# write the angle file
p1.writeangle('angle.x')

# write the dihedral file
p1.writedihedral('dihedral.x')

# write the improper file
p1.writeimproper('improper.x')

