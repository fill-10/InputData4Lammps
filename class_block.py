# class_block.py
# Only linear blocks are available here.
# No dendrimers, no rings.

from class_monomer import monomer
#from class_bond import Bond
#from class_angle import Angle

class block:
	def __init__(self, blockname, monomerobject, repeating, monomerlinkvector):
		self.name = blockname
		self.mon = monomerobject
		self.monrepeat = repeating
		self.Natom = monomerobject.Natom * repeating


		# The monlinkbond is the bond object linking the monomers inside the block. 
		# This is not the bond linking the blocks in the polymer
		self.monlinkvector = monomerlinkvector


		self.deltaV = [\
		self.mon.deltaV[0]*self.monrepeat +self.monlinkvector[0]* (self.monrepeat-1),\
		self.mon.deltaV[1]*self.monrepeat +self.monlinkvector[1]* (self.monrepeat-1),\
		self.mon.deltaV[2]*self.monrepeat +self.monlinkvector[2]* (self.monrepeat-1) ]


		# Set up two attributes to store the pre and suc bonds linking adjacent blocks.
		# This is the inter-block bond
		self.prebond = None
		self.sucbond = None # This sucbond (succeed bond) is not used in linear polymers.

                self.prelinkvector = [1.0,0.0,0.0]



		# Angle information, the two attributes should be the angle objects
		# These are the inter-block angles
		self.headangle_list = []
                self.headdih_list = []


		# check angle
		if self.mon.first_index == self.mon.last_index: # If the monomer's backbone only has one atom
			self.monlinkangle_tail = None # Delete the tail angle constrain because it is reluctant, or, would be contradict to the head angle
			if self.monrepeat <2:
				self.tailangle = None
				# If the block only has one atom in the backbone,
				# i.e. one monomer and one-atom backbone,
				# the inter-block tail angle constrain should be deleted because it would be reluctant.
