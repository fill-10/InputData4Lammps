# bead.py
# bead
class bead:
	def __init__(self, beadname='H', location=[0.,0.,0.],beadtypenumber=0):
		self.name = beadname
		self.coor = location
		self.tpnumber = beadtypenumber
                self.mass = 1.0
