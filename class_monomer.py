# monomer.py
from class_bead import bead
from math import sqrt
#from class_bond import Bond
#from class_angle import Angle
#from class_dihedral import Dihedral
# NOTE: improper constrain is not defined yet.

# calculate the length vector of the monomer, from first pointing to the last bead
def deltaVector(bead2, bead1):
    return [bead2.coor[0] - bead1.coor[0] , \
        bead2.coor[1] - bead1.coor[1] , \
		bead2.coor[2] - bead1.coor[2] ]

class monomer:
    def __init__(self,beadthread = []):
        self.beadlist = beadthread
        self.Natom = len(self.beadlist)
        # first and last bead, by default
        # these can be modified in the main script for a more complicated topology
        self.first_idx = 0
        self.last_idx  = self.Natom - 1
        # Calculate the deltaV vector and length of monomer
        self.deltaV = deltaVector(self.beadlist[self.first_idx], self.beadlist[self.last_idx])
        self.lmon = sqrt( self.deltaV[0]*self.deltaV[0] \
                        + self.deltaV[1]*self.deltaV[1] \
                        + self.deltaV[2]*self.deltaV[2] )
        #
        # Create an empty bond list and an empty anglelist, inside the monomer
        self.bondlist = []
        self.anglelist = []
        self.dihedrallist = []
        self.improperlist = []
        #
        # The input index starts from 0.
    def update_first(self,idx):
        self.first_idx = idx
        self.deltaV = deltaVector( self.beadlist[self.last_idx], \
                                   self.beadlist[self.first_idx] )
        self.lmon = sqrt( self.deltaV[0]*self.deltaV[0] \
                        + self.deltaV[1]*self.deltaV[1] \
                        + self.deltaV[2]*self.deltaV[2] )
    def update_last(self,idx):
        self.last_idx = idx
        self.deltaV = deltaVector( self.beadlist[self.last_idx], \
                                   self.beadlist[self.first_idx] )
        self.lmon = sqrt( self.deltaV[0]*self.deltaV[0] \
                        + self.deltaV[1]*self.deltaV[1] \
                        + self.deltaV[2]*self.deltaV[2] )

