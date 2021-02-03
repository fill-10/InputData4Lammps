# This file is the bead library.
# Define all the beads here.
# Format:
# bead(name, coordinates, typenumber)
#
from class_bead import bead

b1 = bead('S', [0., 1.0, 0.],2)

b2 = bead('B', [0., 0.0, 0.],1)

b1p = bead('S', [0.,-1.0, 0.],2)

b4 = bead('M', [0, -2.0, 0.],4)

b5 = bead('A', [0.,-3.0, 0.],5)


b3 = bead('E',[0.0, 0.0, 0.0], 3)

b6 = bead('O', [0., 0., 0.], 6)

b7 = bead('W', [0., 0., 0.],7)
#
#
PPO = bead('PPO',[0.0, 0., 0.],1)
PPO1 = bead('PPO',[1.0, 0., 0.],1)
PPO2 = bead('PPO',[2.0, 0., 0.],1)
PPO3 = bead('PPO',[3.0, 0., 0.],1)
PA = bead('PA', [1.0, -1., 0.],2)

#
F3 = bead('F3', [0., 0.0, 0.],1)
F31 = bead('F3', [1., 0.0, 0.],1)
F32 = bead('F3', [2., 0.0, 0.],1)
F33 = bead('F3', [3., 0.0, 0.],1)
F34 = bead('F3', [4., 0.0, 0.],1)
FO = bead('FO', [0.,-1.0, 0.],2)
FO1= bead('FO', [0.,-2.0, 0.],2)
SF = bead('SF', [0.,-3.0, 0.],3)

#
BS = bead('Bs', [0., 0.0, 0.],1)
PH = bead('Ph', [0., 1.0, 0.],2)
PHH = bead('Ph',[0.,-1.0, 0.],2)
BS1 = bead('Bs', [1., 0.0, 0.],1)
PH1 = bead('Ph', [1., 1.0, 0.],2)
PH1H = bead('Ph',[1.,-1.0, 0.],2)
BS2 = bead('Bs', [2., 0.0, 0.],1)
PH2 = bead('Ph', [2., 1.0, 0.],2)
PH2H = bead('Ph',[2.,-1.0, 0.],2)
BS3 = bead('Bs', [3., 0.0, 0.],1)
PH3 = bead('Ph', [3., 1.0, 0.],2)
PH3H = bead('Ph',[3.,-1.0, 0.],2)

MIm = bead('Im',[0.,-2.0, 0.],3)
#


W = bead('W', [0., 0., 0.],4)
