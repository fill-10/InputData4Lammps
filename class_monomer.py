# monomer.py
from class_bead import bead
from math import sqrt
from class_bond import Bond
from class_angle import Angle
from class_dihedral import Dihedral
# NOTE: improper constrain is not defined yet.

# calculate the length vector of the monomer, from first pointing to the last bead
def deltaVector(bead2, bead1):
	return [bead2.coor[0] - bead1.coor[0] , \
		bead2.coor[1] - bead1.coor[1] , \
		bead2.coor[2] - bead1.coor[2] ]
def mlength(deltaV):
	return sqrt(pow(deltaV[0],2) + pow(deltaV[1],2) + pow(deltaV[2],2))

class monomer:
	def __init__(self,monomername,monomerindex,beadthread):
		self.name= monomername
		self.index= monomerindex
		
		# input beadlist by hand
		self.beadlist = []
		
		for bd in beadthread:
			# only append the pointers
			self.beadlist.append(bd)
		#
		self.Natom = len(self.beadlist)

		# first and last bead are assigned to the first and last elements in the monomer list, by default
		# these can be modified in the main script, in order to create a more complicated topology
		if len(self.beadlist) > 0 :
			self.first =  self.beadlist[0]
			self.first_index = 1
			#
			# only one bead
			if len(self.beadlist) <2 :
				self.last_index = 1
				self.last = self.beadlist[0]
				self.second_index = 1
				self.second = self.beadlist[0]
				self.lastbut1_index = 1
				self.lastbut1 = self.beadlist[0]
			#
			# 2 beads
			elif len(self.beadlist) <3 :
				self.last_index = 2
				self.last = self.beadlist[1]
				self.second_index = 2
				self.second = self.beadlist[1]
				self.lastbut1_index = 1
				self.lastbut1 = self.beadlist[0]
			#
			# 3 beads
			elif len(self.beadlist) <4 :
				self.last = self.beadlist[2]
				self.last_index = 3
				self.second_index = 2
				self.second = self.beadlist[1]
				self.lastbut1_index = 2
				self.lastbut1 = self.beadlist[1]
			# 4 or more beads
                        else:
				self.last_index = self.Natom
				self.last = self.beadlist[-1]
				self.second_index = 2
				self.second = self.beadlist[self.second_index-1]
				self.lastbut1_index = self.Natom-1
				self.lastbut1 = self.beadlist[self.lastbut1_index-1]
				#
		# Calculate the deltaV vector and length of monomer
		self.deltaV = deltaVector(self.last,self.first)
		self.lmon = mlength(self.deltaV)
		#
		# Create an empty bond list and an empty anglelist, inside the monomer
		self.bondlist = []
		self.anglelist = []
                self.dihedrallist = []
                self.improperlist = []
	#NOTE: The functions below can introduce the index out of dimension. Use these funcitons with care.
        #NOTE: The input index starts from 1.
        #      Index internal starts from 0.
	def update_first(self,index):
		self.first = self.beadlist[index-1]
		self.first_index = index
		self.deltaV = deltaVector(self.last,self.first)
		self.lmon = mlength(self.deltaV)
	#
	def update_last(self,index):
		self.last = self.beadlist[index-1]
		self.last_index = index
		self.deltaV = deltaVector(self.last,self.first)
		self.lmon = mlength(self.deltaV)
	#
	def update_second(self,index):
		self.second = self.beadlist[index-1]
		self.second_index = index
	#
	def update_lastbut1(self,index):
		self.lastbut1 = self.beadlist[index-1]
		self.lastbut1_index = index
	#
	def add_bead(self,b):
		# Note that the add_bead function does not update the last bead
		self.beadlist.append(b)
		self.Natom += 1
	#
	def set_bond_all(self,bondmap):
		self.bondlist = bondmap
	def add_bond(self, onebond):
		self.bondlist.append(onebond)
	#
	def set_angle_all(self, anglemap):
		self.anglelist = anglemap
	def add_angle(self, oneangle):
		self.anglelist.append(oneangle)
        #
        def set_dihedral_all(dihedralmap):
                self.dihedrallist = dihedramap
        def add_angle(self, onedihedral):
                self.dihedrallist.append(onedihedral)


