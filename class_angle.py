# class_angle.py
# angle data structure


# angle class is abandoned in 2.0 version.

class Angle:
	def __init__(self, angletype, beadindex1, beadindexapex, beadindex3):
		self.i1 = beadindex1
		self.iapex = beadindexapex
		self.i3 = beadindex3
		self.angletp = angletype
