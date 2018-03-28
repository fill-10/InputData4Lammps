import os
from class_polymer import polymer
class data:
    def __init__(self,outputfilename):
        self.nsys = 0 # total beads
        #self.nm = 0 # beads per chain, useless
        #self.nmol = 0 # number of chains, useless
        self.nbond = 0 # bond
        self.nangle = 0 # angle
        self.ndihedral = 0 # dihedral
        self.nimproper = 0

        self.all_molecule_list = []
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
        #
        # set up output pointer
        self.output = open(outputfilename, 'w+')
        #

    def load_chain(self, chainxyzsourcefile):
        try:
            chain = open(chainxyzsourcefile, 'r')
            #read nm from chain.xyz
            chain.seek(0)
            chainfirstline = chain.readline()
            number_of_atoms = int(chainfirstline.split()[0])
            #self.chain.close()
            return (number_of_atoms, chain)
        except:
            return (0, None)
    def load_bond(self, bondsourcefile):
        # read No. of bonds per chain from bond file
        try:
        # read all bonds in chain
            bond = open(bondsourcefile, 'r')
            # calculate nbond, bonds per chain
            bond.seek(0)
            bondlist = bond.readlines()
            nbond = len(bondlist)
            return (nbond, bondlist)
        except:
            return (0, None)
    def load_angle(self,anglesourcefile):
        try:
            angle = open(anglesourcefile, 'r')
            angle.seek(0)
            anglelist = angle.readlines()
            nangle =len(anglelist)
            return ( nangle, anglelist)
        except:
            return (0, None)

    def load_dihedral(self, dihedralsourcefile):
        try:
            dihedral = open(dihedralsourcefile, 'r')
            dihedral.seek(0)
            dihedrallist = dihedral.readlines()
            ndihedral = len(dihedrallist)
            return (ndihedral , dihedrallist)
        except:
            return (0, None)

    def load_improper(self, impropersourcefile):
        try:
            improper = open(impropersourcefile, 'r')
            improper.seek(0)
            improperlist = improper.readlines()
            nimproper = len(improperlist)
            return (nimproper, improperlist)
        except:
            return (0, None)


    def read_single_molecule(self, **kwargs):
        chain_filename = None
        bond_filename = None
        angle_filename = None
        dihedral_filename = None
        improper_filename = None
        natom = 0
        for keywords in kwargs.items():
            if keywords[0] == 'chain' or keywords[0] =='molecule':
                chain_filename = keywords[1]
            if keywords[0] == 'bond':
                bond_filename = keywords[1]
            if keywords[0] == 'angle':
                angle_filename = keywords[1]
            if keywords[0] == 'dihedral' :
                dihedral_filename = keywords[1]
            if keywords[0] == 'improper':
                improper_filename = keywords[1]
            if keywords[0] == 'Natom' or keywords[0] == 'natom' or keywords[0] =='number_of_atoms':
                natom = keywords[1]
        return (natom, chain_filename, bond_filename, angle_filename, dihedral_filename, improper_filename)

    def load_all_molecules(self, *some_molecules):
        for molecule in some_molecules:
            self.all_molecule_list.append( self.read_single_molecule(**molecule) )
    
    def update_totalnumber(self):
        for molecule in self.all_molecule_list:
            self.nangle += self.load_angle(molecule[3])[0] * molecule[0]
            self.nbond +=  self.load_bond(molecule[2])[0] * molecule[0]
            self.ndihedral += self.load_dihedral(molecule[4])[0] * molecule[0]
            self.nimproper += self.load_improper(molecule[5])[0] * molecule[0]
            # no need to update self.nsys
            # it's from mixture.xyz






    def load_mixture(self, mixturexyzsourcefile):
        self.mixture = open(mixturexyzsourcefile, 'r')
        # read nsys from mixture.xyz
        self.mixture.seek(0)
        self.mixturefirstline = self.mixture.readline()
        self.nsys = int(self.mixturefirstline.split()[0])



    # useless in this version
    def load_nmol(self, packmolinput):
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

        
    
    # write the head lines of lammps data.
    def writehead(self):
        # some parameters 
        self.output.write('LAMMPS Description \n\n')
        self.output.write(4*' '+ '%8d' %self.nsys + 4*' '+'atoms\n')
        self.output.write(4*' '+ '%8d' %(self.nbond) + 4*' '+'bonds\n')
        self.output.write(4*' '+ '%8d' %(self.nangle) + 4*' '+'angles\n')
        self.output.write(4*' '+ '%8d' %(self.ndihedral) + 4*' '+'dihedrals\n')
        self.output.write(4*' '+ '%8d' %(self.nimproper) + 4*' '+'impropers\n\n')
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
        self.mixture.seek(0)
        self.mixture.readline()
        self.mixture.readline()
        num_prev_mol = 0
        num_prev_atom = 0
        for molecule in self.all_molecule_list:
            natom_per_chain = int( self.load_chain( molecule[1] )[0]   ) # number of atoms per chain
            print(molecule)
            for mol_num_counter in range(0, molecule[0]): # loop over a single type of molecule
                for j in range(0, natom_per_chain ) :
                    atomnamenumber = 0 # claim atom index
                    cl= self.mixture.readline().split() # currentline
                    #
                    # translate:
                    for beadbead in bnlist:
                        if cl[0] == beadbead.name:
                            atomnamenumber = beadbead.tpnumber
                            break
                    #
                    self.output.write('%8d' %(num_prev_atom + mol_num_counter*natom_per_chain +j+1) + ' ' \
                    + '%8d' %(num_prev_mol+mol_num_counter+1) +2*' ' \
                    + '%4d' %atomnamenumber  +' '\
                    +'    0.0000'+' ' \
                    + '%11.4f' %float(cl[1]) + '%11.4f' %float(cl[2]) + '%11.4f' %float(cl[3]) + '\n')
            
            num_prev_mol  += molecule[0]
            num_prev_atom += molecule[0] * natom_per_chain
            try: del atomnamenumber
            except: pass
            #

        self.output.write('\n')

    def writebond(self):
        self.output.write('Bonds\n\n')
        num_prev_bond = 0
        num_prev_atom = 0
        for molecule in self.all_molecule_list:
            nbond_per_chain , current_bond_list= self.load_bond( molecule[2] ) 
            # number of atoms per chain, bond list of curent molecule (raw format)
            natom_per_chain = self.load_chain(molecule[1]  )[0]
            for mol_num_counter in range(0, molecule[0]): # loop over a single type of molecule
                for j in range(0, nbond_per_chain ) :
                    cl=  current_bond_list[j].split() # currentline
                    bondindex = 0 # claim bond index. 0 is an abnormal value.
                    ##################################################################
                    #
                    # The bond type is the FIRST char in each line in the file of bond.x
                    #
                    ##################################################################
                    bondindex = int(cl[0])
                    self.output.write('%8d' %(num_prev_bond+ mol_num_counter * nbond_per_chain +j+1)+' '\
                    + '%4d'%bondindex +' '\
                    +'%7d' %(int(cl[1])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                    +'%7d' %(int(cl[2])+num_prev_atom + mol_num_counter* natom_per_chain)+'\n')
            num_prev_atom += molecule[0] * natom_per_chain
            num_prev_bond += molecule[0] * nbond_per_chain

        self.output.write('\n')
        print("Bonds written")


    def writeangle(self):
        self.output.write('Angles\n\n')
        num_prev_angle = 0
        num_prev_atom = 0
        for molecule in self.all_molecule_list:
            nangle_per_chain , current_angle_list= self.load_angle( molecule[3] ) 
            # number of atoms per chain, angle list of curent molecule (raw format)
            natom_per_chain = self.load_chain(molecule[1]  )[0]
    
            for mol_num_counter in range(0, molecule[0]): # loop over a single type of molecule
                for j in range(0, nangle_per_chain) :
                    cl=  current_angle_list[j].split() # currentline
                    angleindex = 0 # claim bond index. 0 is an abnormal value.
            ##################################################################
            #
            # The angle type is the FIRST char in each line in the file of bond.x
            #
            ##################################################################
                    angleindex = int(cl[0])
                    self.output.write('%8d' %(num_prev_angle+ mol_num_counter * nangle_per_chain +j+1)+' '\
                    + '%4d'%angleindex +' '\
                    +'%7d' %(int(cl[1])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                    +'%7d' %(int(cl[2])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                    +'%7d' %(int(cl[3])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                   +'\n')
            num_prev_atom += molecule[0] * natom_per_chain
            num_prev_angle += molecule[0] * nangle_per_chain

        self.output.write('\n')
        print("Angles written")

    def writedihedral(self):
        self.output.write('Dihedrals\n\n')
        num_prev_dihedral = 0
        num_prev_atom = 0
        for molecule in self.all_molecule_list:
            ndihedral_per_chain , current_dihedral_list= self.load_dihedral( molecule[4] ) 
            # number of atoms per chain, dihedral list of curent molecule (raw format)
            natom_per_chain = self.load_chain(molecule[1]  )[0]
            for mol_num_counter in range(0, molecule[0]): # loop over a single type of molecule
                for j in range(0, ndihedral_per_chain) :
                    cl=  current_dihedral_list[j].split() # currentline
                    dihedralindex = 0 # claim index. 0 is an abnormal value.
                    dihedralindex = int(cl[0])
                    self.output.write('%8d' %(num_prev_dihedral+ mol_num_counter * ndihedral_per_chain +j+1)+' '\
                    + '%4d'%dihedralindex +' '\
                    +'%7d' %(int(cl[1])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                    +'%7d' %(int(cl[2])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                    +'%7d' %(int(cl[3])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                    +'%7d' %(int(cl[4])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                    +'\n')
            num_prev_atom += molecule[0] * natom_per_chain
            num_prev_dihedral += molecule[0] * ndihedral_per_chain

        self.output.write('\n')
        print("Dihedrals written")


    def writeimproper(self):
        self.output.write('Impropers\n\n')
        num_prev_improper = 0
        num_prev_atom = 0
        for molecule in self.all_molecule_list:
            nimproper_per_chain , current_improper_list= self.load_improper( molecule[5] ) 
            # number of atoms per chain, improper list of curent molecule (raw format)
            natom_per_chain = self.load_chain(molecule[1]  )[0]
            for mol_num_counter in range(0, molecule[0]): # loop over a single type of molecule
                for j in range(0, nimproper_per_chain) :
                    cl=  current_improper_list[j].split() # currentline
                    improperindex = 0 # claim index. 0 is an abnormal value.
                    improperindex = int(cl[0])
                    self.output.write('%8d' %(num_prev_improper+ mol_num_counter * nimproper_per_chain +j+1)+' '\
                    + '%4d'%improperindex +' '\
                    +'%7d' %(int(cl[1])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                    +'%7d' %(int(cl[2])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                    +'%7d' %(int(cl[3])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                    +'%7d' %(int(cl[3])+num_prev_atom + mol_num_counter* natom_per_chain)+' '\
                    +'\n')
            num_prev_atom += molecule[0] * natom_per_chain
            num_prev_improper += molecule[0] * nimproper_per_chain

        self.output.write('\n')
        print("Impropers written")


