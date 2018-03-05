# class_packmolgen.py
# generate packmol input
import os
import subprocess


class packmolgen:
	def __init__(self, filename, boxsize, components, filetype='xyz', tolerance= 0.1):
		self.xlo = 0.
		self.ylo = 0.
		self.zlo = 0.
		self.xhi = boxsize[0]
		self.yhi = boxsize[1]
		self.zhi = boxsize[2]
		self.tol = tolerance
		self.ftype = filetype
		self.allcomp = components
		self.inpfilename = filename
	
	def packmolinpgen(self, outputmixturefilename):
		self.packmolf = open(self.inpfilename, 'w+')
		self.packmolf.write('tolerance  ' + str(self.tol) + '\n\n')
		self.packmolf.write('filetype  ' + self.ftype + '\n\n')
		self.packmolf.write('output  ' + outputmixturefilename + '\n\n')
		for comp in self.allcomp:
			self.packmolcomponent(comp[0], comp[1])
		self.packmolf.close()

	def packmolcomponent(self, compfilename, number):
		self.packmolf.write('structure  ' + compfilename +'\n')
		self.packmolf.write('  number  ' + str(number) + '\n')
		self.packmolf.write('  inside box  ' +\
			'%10.4f' %self.xlo + '%10.4f' %self.ylo +'%10.4f' %self.zlo +2*' ' +\
			'%10.4f' %self.xhi + '%10.4f' %self.yhi + '%10.4f' %self.zhi +'\n')
		self.packmolf.write('end structure \n\n')
	
	def runpackmol(self, packmolcommand):
		fullcommand = packmolcommand+'<'+self.inpfilename
		subprocess.call(fullcommand)
