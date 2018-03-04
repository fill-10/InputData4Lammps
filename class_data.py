import os
from class_polymer import polymer
class data:
	def __init__(self,outputfilename,mixturexyzsourcefile=None, chainxyzsourcefile=None, bondsourcefile=None, anglesourcefile = None, packmolinput = None):
		
		self.nsys = 0 # total beads
		self.nm = 0 # beads per chain
		self.nmol = 0 # number of chains
		self.nbond = 0 # bond per chain
		self.nangle = 0 # angle per chain
		self.ndihedral = 0
		self.nimproper = 0
		
		# set these values manually. No. of types 
		self.tatom = 3
		self.tbond = 5
		self.tangle = 5
		self.tdihedral = 0
		self.timproper = 0

		# box dimensions
		self.xlo = 0.
		self.ylo = 0.
		self.zlo = 0.
		self.xhi = 120.
		self.yhi = 120.
		self.zhi = 120.
		
		# read total No. of beads from mixture.xyz. it is in its first line.
		if mixturexyzsourcefile:
			self.mixture = open(mixturexyzsourcefile, 'r')
			# read nsys from mixture.xyz
			self.mixture.seek(0)
			self.mixturefirstline = self.mixture.readline()
			self.nsys = int(self.mixturefirstline.split()[0])
		else: self.mixture = None
		
		# read No. of beads per chain from chain.xyz
		if chainxyzsourcefile:
			self.chain = open(chainxyzsourcefile, 'r')
			#read nm from chain.xyz
			self.chain.seek(0)
			self.chainfirstline = self.chain.readline()
			self.nm = int(self.chainfirstline.split()[0])
			#self.chain.close()
		else: self.chain = None
		
		# read No. of chains from packmol input
		if packmolinput:
			self.packmol = open(packmolinput,'r')
			packmoldata = self.packmol.readlines()
			# read the number of chains
			for lineindex in range(0, len(packmoldata)):
				if packmoldata[lineindex].find('chain') > -1 :
					self.nmol = int(packmoldata[lineindex+1].split()[1])
					break
			# delete this big local variable. thanks to the small-sized packmol input, or this method would suck...
			#del packmoldata
			#del lineindex
		else:
			self.packmol = None

		# read No. of bonds per chain from bond file
		if bondsourcefile:
			# read all bonds in chain
			self.bond = open(bondsourcefile, 'r')
			# calculate nbond, bonds per chain
			self.bond.seek(0)
			self.bondlist = self.bond.readlines()
			self.nbond = len(self.bondlist)
		else: self.bond = None
		
		# read No. of angles per chain from angle file
		if anglesourcefile:
			self.angle = open(anglesourcefile, 'r')
			self.angle.seek(0)
			self.anglelist = self.angle.readlines()
			self.nangle =len(self.anglelist)
		else: self.angle = None

		#
		# set up output pointer
		self.output = open(outputfilename, 'w+')
		#

	# write the head lines of lammps data.
	def writehead(self):
		# some parameters 
		self.output.write('LAMMPS Description \n\n')
		self.output.write(4*' '+ '%8d' %self.nsys + 4*' '+'atoms\n')
		self.output.write(4*' '+ '%8d' %(self.nbond*self.nmol) + 4*' '+'bonds\n')
		self.output.write(4*' '+ '%8d' %(self.nangle*self.nmol) + 4*' '+'angles\n')
		self.output.write(4*' '+ '%8d' %(self.ndihedral*self.nmol) + 4*' '+'dihedrals\n')
		self.output.write(4*' '+ '%8d' %(self.nimproper*self.nmol) + 4*' '+'impropers\n\n')
		
		self.output.write(4*' '+ '%8d' %self.tatom + 4*' '+'atom types\n')
		self.output.write(4*' '+ '%8d' %self.tbond + 4*' '+'bond types\n')
		self.output.write(4*' '+ '%8d' %self.tangle + 4*' '+'angle types\n')
		self.output.write(4*' '+ '%8d' %self.tdihedral + 4*' '+'dihedral types\n')
		self.output.write(4*' '+ '%8d' %self.timproper + 4*' '+'improper types\n\n')
		
		# box size
		self.output.write('%21.14f' %self.xlo + 5*' '+'%21.14f' %self.xhi+6*' ' + 'xlo'+' '+'xhi\n')
		self.output.write('%21.14f' %self.ylo + 5*' '+'%21.14f' %self.yhi+6*' ' + 'ylo'+' '+'yhi\n')
		self.output.write('%21.14f' %self.zlo + 5*' '+'%21.14f' %self.zhi+6*' ' + 'zlo'+' '+'zhi\n\n')

		# Masses
		self.output.write('Masses \n\n')
		for i in range(0, self.tatom):
			self.output.write(4*' '+ '%8d' %(i+1) + 2*' '+ '%11.8f' %1.+'\n')
		self.output.write('\n')
	
	def writeatom(self, bnlist):
		
		self.output.write('Atoms\n\n')

		# need to skip the first 2 lines in the mixture.xyz
		# but in the __init__, 1 line has already been read
		# thus, here only need to skip 1 line
		self.mixture.readline()
		
		for i in range(0, self.nmol):
			for j in range(0, self.nm):
				atomnamenumber = 0 # claim atom index
				cl= self.mixture.readline().split() # currentline
				for beadbead in bnlist:
					if cl[0] == beadbead.name:
						atomnamenumber = beadbead.namenumber
						break
				self.output.write('%8d' %(i*self.nm+j+1) + ' ' + '%8d' %(i+1) +2*' ' + '%4d' %atomnamenumber \
				+' '+'    0.0000'+' ' + '%11.4f' %float(cl[1]) + '%11.4f' %float(cl[2]) + '%11.4f' %float(cl[3]) + '\n')
		# then, need to write the water and OH beads~~~ 
		for i in range(self.nmol*self.nm, self.nsys):
			atomnamenumber = 0
			cl = self.mixture.readline().split()
			for beadbead in bnlist:
				if cl[0] == beadbead.name:
					atomnamenumber = beadbead.namenumber
					break
			self.output.write('%8d' %(i+1) + ' ' + '%8d' %(self.nmol + i-self.nmol*self.nm+1) +2*' ' + '%4d' %atomnamenumber \
			+' '+'    0.0000'+' ' + '%11.4f' %float(cl[1]) + '%11.4f' %float(cl[2]) + '%11.4f' %float(cl[3]) + '\n')
		self.output.write('\n')

	def writebond(self):
		#
		# the bondnamelist is useless in the new version
		#

		self.output.write('Bonds\n\n')
		for i in range (0, self.nmol):
			for j in range(0, self.nbond):
				cl = self.bondlist[j].split()

				bondindex = 0 # claim bond index. 0 is an abnormal value.

				##################################################################
				#
				# This section uses bead name to recognize the bond type.
				# This method is discarded in this case.
				# The bond type is the last char in each line in the new version of bond.x
				#
				##################################################################

				##################################################################
				# 
				# New version of bond type recognition.
				# It is simpler.
				bondindex = int(cl[4])

				self.output.write('%8d' %(i*self.nbond+j+1)+' '+ '%4d'%bondindex \
				+' '+'%7d' %(int(cl[0])+i*self.nm)+' '+'%7d' %(int(cl[1])+i*self.nm)+'\n')

		self.output.write('\n')


	def writeangle(self):
		self.output.write('Angles\n\n')
		for i in range (0, self.nmol):
			for j in range(0, self.nangle):
				cl = self.anglelist[j].split()
				angleindex = 0 # claim angle index. 0 is an abnormal value.


				############################################################
				#
				# Same as bond, the old version is based on the names of bead.
				# The new version directly reads the angle type from angle.x.
				#
				###############################################################

				angleindex = int(cl[6])

				self.output.write('%8d' %(i*self.nangle+j+1)+' '+ '%4d'%angleindex \
				+' ' +'%7d' %(int(cl[0])+i*self.nm)+' '+'%7d' %(int(cl[1])+i*self.nm)+' '+'%7d' %(int(cl[2])+i*self.nm)+'\n')
        def writedihedral(self):
            pass
