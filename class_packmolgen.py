# class_packmolgen.py
# generate packmol input
import os
import subprocess


class packmolgen:
	def __init__(self, **kwargs):
	        self.xlo = 0.
		self.ylo = 0.
		self.zlo = 0.
		self.xhi = 100.
		self.yhi = 100.
		self.zhi = 100.
		self.tol = 0.1
		self.ftype = 'xyz'
		self.allcomp = []
		self.inpfilename = 'packmol.inp'
                for kw in kwargs.items():
                    if kw[0] == 'boxsize':
                        self.xhi = kw[1][0]
                        self.yhi = kw[1][1]
                        self.zhi = kw[1][2]
                    if kw[0] == 'molecules':
                        self.molconfig = kw[1]
                        for molecule in self.molconfig:
                            self.allcomp.append([molecule['chain'], molecule['Natom']])
                    if kw[0] == 'tolerance':
                        self.tol = kw[1]
                    if kw[0] == 'filetype':
                        self.ftype = kw[1]
                    if kw[0] == 'inputfilename':
                        self.inpfilename = kw[1]

	
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
