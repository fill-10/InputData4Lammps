# class_bond.py
# bond data structure
# Example: bond1 = Bond(1, 2, 1, [0., 1., 0.])

class Bond:
	def __init__(self, bondtypeindex, beadindex1, beadindex2, bondvector=[1.0, 0.0, 0.0]):
		self.i1 = beadindex1
		self.i2 = beadindex2
		if beadindex1 > beadindex2:
			self.i1 = beadindex2
			self.i2 = beadindex1
		self.btype = bondtypeindex
		self.vector = bondvector
