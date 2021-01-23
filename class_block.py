# class_block.py
# Only linear blocks are available here.
# No dendrimers, no rings.

from class_monomer import monomer

class block:
    def __init__(self, monomerobj, repeat, monlinkvector):
        self.mon = monomerobj
        self.monrep = repeat
        self.Natom = self.mon.Natom * self.monrep

        # self.monlink is mon-> mon 
        # It is not the bond linking the blocks.
        self.monlink = monlinkvector

        self.deltaV = [\
        self.mon.deltaV[0]*self.monrep +self.monlink[0]* (self.monrep-1), \
        self.mon.deltaV[1]*self.monrep +self.monlink[1]* (self.monrep-1), \
        self.mon.deltaV[2]*self.monrep +self.monlink[2]* (self.monrep-1) ]

        # Set up two attributes to store the pre and suc bonds linking adjacent blocks.
        # This is the inter-block bond
        self.prebond = []
        self.prelinkvector = [0.0,0.0,0.0]
        # Angle information, the two attributes should be the angle objects
        # These are the inter-block angles
        self.headangle_list = []
        self.headdih_list = []
        self.headimproper_list = []
