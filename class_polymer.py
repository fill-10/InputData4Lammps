# class_polymer.py
#
#
from class_block import block
#from class_bond import Bond

import os


class polymer:
	def __init__(self, polymername, blocklist):
		# The order of the blocks matters. They are linked one after one in the order of the input.
		self.name = polymername
		# create an empty list of blocks
		self.block_list = []
		self.TotalNatom = 0
		# add blocks to the list
		for i in range(0,len(blocklist)):
			current_block = blocklist[i]
			if i == 0:
				current_block.prebond = None
			if i == len(blocklist)-1:
				current_block.sucbond = None
			self.block_list += [current_block]
			# calculate the total number of beads in the polymer
			self.TotalNatom += current_block.Natom

	#
	#################################################################
	# This function, writexyz, writes the atoms/beads into file.    #
	# 
	def writexyz(self,filename='chain1.xyz'):
		# creat/overwrite the xyz in read and write mode
		self.f= open(filename, 'w+')

		self.f.write(str(self.TotalNatom)+'\n\n')
		#
		# Use a variable pre_block_Vector to count for the size vector of all previous blocks
		# NOTE: Don't directly use pre_block_Vector = self.block_list[0].mon.first.coordinates!
		# In Python, the operator '=' passes the reference!
		# So, you need to claim the pre_block_Vector=[], separately.
		# Then append the self.block_list[0].mon.first.coordinates to this new vector!
		# Or, it would share the same mem space with self.block_list[0].mon.first.coordinates,
		# and the self.block_list[0].mon.first.coordinates would be modified when accumulate pre_block_Vector!
		pre_block_Vector = []
		pre_block_Vector += self.block_list[0].mon.first.coor

		#
		for i in range(0, len(self.block_list)):
			if i >0:
				pre_block_Vector[0] += self.block_list[i-1].deltaV[0] + self.block_list[i].prelinkvector[0]
				pre_block_Vector[1] += self.block_list[i-1].deltaV[1] + self.block_list[i].prelinkvector[1]
				pre_block_Vector[2] += self.block_list[i-1].deltaV[2] + self.block_list[i].prelinkvector[2]

			for r in range(0, self.block_list[i].monrepeat):
				for bd in self.block_list[i].mon.beadlist:
					self.f.write\
					( '%s' %bd.name + 3*' ' \
					+ '%6.2f' %( pre_block_Vector[0] + bd.coor[0] + r * ( self.block_list[i].mon.deltaV[0] + self.block_list[i].monlinkvector[0] )  ) + 3*' ' \
					+ '%6.2f' %( pre_block_Vector[1] + bd.coor[1] + r * ( self.block_list[i].mon.deltaV[1] + self.block_list[i].monlinkvector[1] ) ) + 3*' ' \
					+ '%6.2f' %( pre_block_Vector[2] + bd.coor[2] + r * ( self.block_list[i].mon.deltaV[2] + self.block_list[i].monlinkvector[2] ) ) + 3*' ' \
					+ '\n' \
					)

		self.f.close()
	#
	#
	########################################################
	# This function writes the bond information into file. #
	#
	def writebond(self,filename):
		self.f = open(filename, 'w+')
		
		preNatom= 0 #set up a counter to count the number of atoms in all previously covered blocks.
		
		
		for i in range(0, len(self.block_list)):

			cbk = self.block_list[i] # current block

			# bond linking the previous blocks, starting from the second block:
			if i>0:
				pbk = self.block_list[i-1] # read the previous block
				if cbk.prebond:
					self.f.write(\
                                        '%8d' %cbk.prebond[0] \
					+ '%12d' %( (preNatom - pbk.Natom) + (pbk.monrepeat-1) * pbk.mon.Natom  + pbk.mon.last_index) \
					+ '%12d' %(  preNatom                                                   + cbk.mon.first_index) \
					+ '%8s' %pbk.mon.last.name \
					+ '%8s' %cbk.mon.first.name \
                                        + '\n')


			for r in range(0, cbk.monrepeat):
				# bonds inside one monomer (normal)
				for j in range(0, len(cbk.mon.bondlist)):
                                    # regular bond inside monomer
                                    if  (None not in cbk.mon.bondlist[j]) \
                                    and not any(bn <= 0 for bn in cbk.mon.bondlist[j]) :
					self.f.write(\
					'%8d' %cbk.mon.bondlist[j][0] \
					+ '%12d'   %(preNatom + r * cbk.mon.Natom + cbk.mon.bondlist[j][1]) \
					+ '%12d' %(preNatom + r * cbk.mon.Natom + cbk.mon.bondlist[j][2]) \
					+ '%8s' %cbk.mon.beadlist[cbk.mon.bondlist[j][1]-1].name \
					+ '%8s' %cbk.mon.beadlist[cbk.mon.bondlist[j][2]-1].name \
					+ '\n')
				# bond connecting monomers
                                # There is None or 0 in bond
                                    elif r>0:
                                        if not cbk.mon.bondlist[j][1]:
                                            tmp1beadname = cbk.mon.last.name
                                            tmp1beadindex = cbk.mon.last_index
                                        else:
                                            tmp1beadindex = -cbk.mon.bondlist[j][1]
                                            tmp1beadname = cbk.mon.beadlist[tmp1beadindex-1].name
                                        if not cbk.mon.bondlist[j][2]:
                                            tmp2beadname = cbk.mon.first.name
                                            tmp2beadindex = cbk.mon.first_index
                                        else:
                                            tmp2beadindex = -cbk.mon.bondlist[j][2]
                                            tmp2beadname = cbk.mon.beadlist[-tmp2beadindex-1].name
                                        
					self.f.write(\
					'%8d' %cbk.mon.bondlist[j][0] \
                                        +'%12d' %(  preNatom + (r-1) * cbk.mon.Natom + tmp1beadindex ) \
					+ '%12d' %(preNatom +  r    * cbk.mon.Natom + tmp2beadindex) \
					+ '%8s' %tmp1beadname\
					+ '%8s' %tmp2beadname \
					+  '\n' )

			preNatom += cbk.Natom

		self.f.close()
	#
	#
	#
	def writeangle(self,filename):
            self.f=open(filename, 'w+')
            preNatom= 0 #set up a counter to count the number of atoms in all previously covered blocks.
            for i in range(0, len(self.block_list)):
                cbk = self.block_list[i] # current block
		################################################
		##
		## Write the angle linking the previous block.  #
		##
		################################################
                
                # translate:

		if i>0:
                    pbk = self.block_list[i-1] # read the previous block
                    if len(cbk.headangle_list)>0:
                        for m in range(0, len(cbk.headangle_list) ):
                            #
                            # if the apex is in cbk
                            if cbk.headangle_list[m][2]>0 :
                                ## some translation:
                                if not cbk.headangle_list[m][1]:
                                    cbk.headangle_list[m][1] = -pbk.mon.last_index
                                #
                                # trouble maker bead name... patched.
                                tmp3bdnm =  cbk.mon.beadlist[cbk.headangle_list[m][3]-1].name
                                
                                # if only a single bead in the backbone of cbk monomer
                                if cbk.mon.lmon ==0.0:
                                    #
                                    # if cbk has more than 1 monomer
                                    if cbk.monrepeat >1 :
                                    #
                                        if cbk.headangle_list[m][2] == cbk.headangle_list[m][3]:
                                            cbk.headangle_list[m][3] += cbk.mon.Natom
                                        elif (not cbk.headangle_list[m][3]):
                                            cbk.headangle_list[m][3] = cbk.mon.first_index + cbk.mon.Natom
                                        else:
                                            print("illegal angle constrain in")
                                            print(self.block_list[i])
                                            return -1

                                    #
                                    if (cbk.monrepeat ==1) :
                                    # The very special case: bk1-bk2-bk3.
                                    # Now doing bk2, bk2 only has one monomer and the monomer has the single-atom backbone.
                                    # The bk2 has the headangle constrain. i.e. the only backbone atom is the apex of the angle.
                                    # and there is a following block
                                        if (i< len(self.block_list)-1): # if not the last block.
                                            cbk.headangle_list[m][3] += cbk.Natom
                                            tmp3bdnm = self.blocklist[i+1].mon.first.name
                                        else:
                                            print("illegal angle constrain in")
                                            print(i)
                                            return -1

                                # write
                                self.f.write(\
                                '%8d'    %cbk.headangle_list[m][0] \
				+'%12d'  %( (preNatom - pbk.Natom) + (pbk.monrepeat-1) * pbk.mon.Natom  - cbk.headangle_list[m][1] ) \
		                + '%12d' %(  preNatom                                                   + cbk.headangle_list[m][2] ) \
			        + '%12d' %(  preNatom                                                   + cbk.headangle_list[m][3] ) \
				+ '%8s'  %pbk.mon.beadlist[-cbk.headangle_list[m][1]-1].name \
				+ '%8s'  %cbk.mon.beadlist[ cbk.headangle_list[m][2]-1].name \
				+ '%8s'  %tmp3bdnm \
                                + '\n' )
                                
                                del tmp3bdnm
                            #
                            # if apex bead in pbk:
                            elif pbk.monrepeat >1: # the 3rd bead is None, 0, negative index
                                # some translation
                                if not cbk.headangle_list[m][2]:
                                    cbk.headangle_list[m][2] = -pbk.mon.last_index
                                if not cbk.headangle_list[m][1]:
                                    cbk.headangle_list[m][1] = -pbk.mon.lastbut1_index
                                #
                                # trouble maker bead name... patched.
                                tmp1bdnm =  pbk.mon.beadlist[-cbk.headangle_list[m][1]-1].name
                                #
                                if cbk.headangle_list[m][1] == cbk.headangle_list[m][2] and \
                                pbk.mon.lmon == 0.0 :
                                    cbk.headangle_list[m][1] += pbk.mon.Natom
                                    if pbk.monrepeat <=1:  #single atom in pbk.mon and no repeat
                                        #
                                        # need to check p2bk...
                                        print( 'illiegal head angular constrain in')
                                        print( i )
                                        return -1
                                # write
                                self.f.write(\
                                '%8d' %cbk.headangle_list[m][0] \
				+'%12d'  %( (preNatom - pbk.Natom) + (pbk.monrepeat-1) * pbk.mon.Natom  - cbk.headangle_list[m][1] ) \
		                + '%12d' %( (preNatom - pbk.Natom) + (pbk.monrepeat-1) * pbk.mon.Natom  - cbk.headangle_list[m][2] ) \
			        + '%12d' %(  preNatom                                                   + cbk.headangle_list[m][3] ) \
				+ '%8s'  %tmp1bdnm \
				+ '%8s'  %pbk.mon.beadlist[-cbk.headangle_list[m][2]-1].name \
				+ '%8s'  %cbk.mon.beadlist[ cbk.headangle_list[m][3]-1].name \
                                + '\n' )

                                del tmp1bdnm


		#################################################
		##
		## Write the angles inside each block.          #
		##
		#################################################
                #
                # This IF is a must.
                # There is the case:
                # only block.headangle_list has constrains,
                # but no angle inside each block,
                # then the headdih will pollute the following lines
                # and cause error.
                if len(cbk.mon.anglelist) == 0:
                    preNatom += cbk.Natom
                    continue
                #
                #
                ###


                # translate :
                for j in range(0, len(cbk.mon.anglelist) ):
                    if cbk.mon.anglelist[j][2] ==0 or cbk.mon.anglelist[j][2] == None:
                        cbk.mon.anglelist[j][2] = cbk.mon.last_index
                        if cbk.mon.anglelist[j][1] ==0 or cbk.mon.anglelist[j][1] ==None:
                            cbk.mon.anglelist[j][1] =  - cbk.mon.lastbut1_index
                    elif cbk.mon.anglelist[j][1] ==0 or cbk.mon.anglelist[j][1] ==None:
                            cbk.mon.anglelist[j][1] =  - cbk.mon.last_index
                    if cbk.mon.anglelist[j][3] ==0 or cbk.mon.anglelist[j][3] == None:
                        cbk.mon.anglelist[j][3] = cbk.mon.first_index
                        # if cbk.mon.anglelist[j][3] != cbk.mon.anglelist[j][1]:
                        #     print "illegal angle constrain in"
                        #     return -1

                del j
                ## translate done
                #

                for r in range(0, cbk.monrepeat):
                # From the second monomer, need to include the angles linking the current and previous monomers.
                # I.e., r from 1 to self.repeat-1 (incld)
                    for j in range(0, len(cbk.mon.anglelist)):
                        try:
                            tmp1bdindex = cbk.mon.anglelist[j][1] 
                            tmp1bdnm = cbk.mon.beadlist[tmp1bdindex-1].name
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp2bdindex = cbk.mon.anglelist[j][2] 
                            tmp2bdnm = cbk.mon.beadlist[tmp2bdindex-1].name
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp3bdindex = cbk.mon.anglelist[j][3] 
                            tmp3bdnm = cbk.mon.beadlist[tmp3bdindex-1].name 
                        except:
                            pass
                        finally:
                            pass
                        #
                        # if the bond is a monomer linking bond,  update:
                        if any(bn <= 0 for bn in cbk.mon.anglelist[j]):
                            if cbk.mon.anglelist[j][2] >0 and cbk.mon.anglelist[j][3]>0:
                                tmp1bdindex = -cbk.mon.anglelist[j][1] - cbk.mon.Natom
                                tmp1bdnm    = cbk.mon.beadlist[ -cbk.mon.anglelist[j][1] -1 ].name
                                if cbk.mon.anglelist[j][3] == cbk.mon.anglelist[j][2]:
                                    tmp3bdindex = cbk.mon.anglelist[j][3] + cbk.mon.Natom
                                elif cbk.mon.lmon==0: # additional check
                                    tmp3bdindex = -cbk.mon.anglelist[j][3] +  ckb.mon.Natom
                                    #tmp3bdnm    = ckb.mon.beadlist[ -cbk.mon.anglelist[j][3] -1]
                                else: 
                                    print("illegal angle constrain:")
                                    print(" monomer does not have a single atom backbon. ")
                                    print("intermonomer angle needs to be modified.")
                                    return -1
                            elif cbk.mon.anglelist[j][2]<0 and cbk.mon.anglelist[j][1]<0:
                                tmp1bdindex  = -cbk.mon.anglelist[j][1] -2*cbk.mon.Natom
                                tmp1bdnm     =  cbk.mon.beadlist[-cbk.mon.anglelist[j][1] -1 ].name
                                tmp2bdindex  = -cbk.mon.anglelist[j][2] -  cbk.mon.Natom
                                tmp2bdnm     =  cbk.mon.beadlist[-cbk.mon.anglelist[j][2] -1 ].name
                            elif cbk.mon.anglelist[j][3] <0:
                                tmp3bdindex  = -cbk.mon.anglelist[j][3] - cbk.mon.Natom
                                tmp3bdnm  =  cbk.mon.beadlist[-cbk.mon.anglelist[j][3] -1 ].name
                            # 2nd round translate done.

                        #write:
                        if r==0 and any(bn <= 0 for bn in cbk.mon.anglelist[j]):
                            # if it is the monomer linking angle
                            # and this is the first monomer
                            # go to the next angle, i.e. j++
                            continue
                        
                        elif r ==1:
                            # second monomer
                            if cbk.mon.anglelist[j][2]<0 and cbk.mon.anglelist[j][1] <0\
                            and cbk.mon.lmon ==0.0:
                                continue
                        elif  r>= cbk.monrepeat-1:
                            # last monomer
                            if cbk.mon.anglelist[j][3] == cbk.mon.anglelist[j][2] >0\
                            and cbk.mon.lmon == 0.0 :
                                continue
                        #print(tmp1bdindex, tmp2bdindex, tmp3bdindex, tmp1bdnm, tmp2bdnm,tmp3bdnm)
                        self.f.write(\
                        '%8d'  %cbk.mon.anglelist[j][0] \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp1bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp2bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp3bdindex) \
                        + '%8s'  %tmp1bdnm \
                        + '%8s'  %tmp2bdnm \
                        + '%8s'  %tmp3bdnm \
                        + '\n')
                #                
                preNatom += cbk.Natom
		
            self.f.close()
        
        def writedihedral(self, filename):
            self.f= open(filename, 'w+')
            preNatom= 0 #set up a counter to count the number of atoms in all previously covered blocks.
            for i in range(0, len(self.block_list)):
                print i
                ###--------------------------###
                ### block linking Dihedral   ###
                ###--------------------------###
                cbk = self.block_list[i]
                if i > 0:
                    pbk = self.block_list[i-1]
                    if len(cbk.headdih_list)>0:
                        for m in range(0, len(cbk.headdih_list) ):
                            tmp1bdindex = cbk.headdih_list[m][1]
                            tmp2bdindex = cbk.headdih_list[m][2]
                            tmp3bdindex = cbk.headdih_list[m][3]
                            tmp4bdindex = cbk.headdih_list[m][4]
                            # translate
                            if cbk.headdih_list[m][3]<0:
                                if pbk.mon.lmon == 0.0:
                                    if pbk.monrepeat >2:
                                        tmp1bdindex = -cbk.headdih_list[m][1] -3*pbk.mon.Natom
                                        tmp2bdindex = -cbk.headdih_list[m][2] -2*pbk.mon.Natom
                                        tmp3bdindex = -cbk.headdih_list[m][3] -  pbk.mon.Natom
                                    else:
                                        print("not supported, sorry.")
                                elif pbk.mon.lastbut1_index == pbk.mon.first_index: # 2 beads in prev blk mon
                                    if pbk.monrepeat >1:
                                        tmp1bdindex = -cbk.headdih_list[m][1] -2*pbk.mon.Natom
                                        tmp2bdindex = -cbk.headdih_list[m][2] -  pbk.mon.Natom
                                        tmp3bdindex = -cbk.headdih_list[m][3] -  pbk.mon.Natom
                                    else:
                                        print("not supported, sorry.")
                            elif cbk.headdih_list[m][2]<0:
                                if pbk.mon.lmon == 0.0:
                                    if pbk.monrepeat >1:
                                        tmp1bdindex = -cbk.headdih_list[m][1] -2*pbk.mon.Natom
                                        tmp2bdindex = -cbk.headdih_list[m][2] -1*pbk.mon.Natom
                                    else:
                                        print("not supported, sorry.")
                                elif pbk.mon.lastbut1_index == pbk.mon.first_index: # 2 beads in prev blk mon
                                    tmp1bdindex = -cbk.headdih_list[m][1] -  pbk.mon.Natom
                                    tmp2bdindex = -cbk.headdih_list[m][2] -  pbk.mon.Natom

                            elif cbk.headdih_list[m][1]<0:
                                tmp1bdindex = -cbk.headdih_list[m][1] - pbk.mon.Natom

                            else:
                                print("illegal head dihedral")
                                print("block:")
                                print(i)
                                return -1
                            # write
                            self.f.write(\
                            '%8d'  %cbk.headdih_list[m][0] \
                            + '%12d' %( preNatom + tmp1bdindex) \
                            + '%12d' %( preNatom + tmp2bdindex) \
                            + '%12d' %( preNatom + tmp3bdindex) \
                            + '%12d' %( preNatom + tmp4bdindex) \
                            + '\n')
                    try:
                        del tmp1bdindex
                    except:
                        pass
                    try:
                        del tmp2dindex
                    except:
                        pass
                    try:
                        del tmp3bdindex
                    except:
                        pass
                    try:
                        del tmp1bdindex
                    except:
                        pass
                ###--------------------------###
                ### dihedral in each block   ###
                ###--------------------------###
                ###
                # This IF is a must.
                # There is the case:
                # only block.headdih_list has dihedra,
                # but no dihedral inside each block,
                # then the headdih will pollute the following lines
                # and cause error.
                if len(cbk.mon.dihedrallist) == 0:
                    preNatom += cbk.Natom
                    continue
                #
                #
                ###
                for r in range(0, cbk.monrepeat):
                    for j in range(0 , len(cbk.mon.dihedrallist)):
                        try:
                            tmp1bdindex = cbk.mon.dihedrallist[j][1] 
                            tmp1bdnm = cbk.mon.beadlist[tmp1bdindex-1].name
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp2bdindex = cbk.mon.dihedrallist[j][2] 
                            tmp2bdnm = cbk.mon.beadlist[tmp2bdindex-1].name
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp3bdindex = cbk.mon.dihedrallist[j][3] 
                            tmp3bdnm = cbk.mon.beadlist[tmp3bdindex-1].name 
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp4bdindex = cbk.mon.dihedrallist[j][4] 
                            tmp4bdnm = cbk.mon.beadlist[tmp3bdindex-1].name 
                        except:
                            pass
                        finally:
                            pass
                        if any(bn <= 0 for bn in cbk.mon.dihedrallist[j]):
                            #
                            # -2 2 3 4 
                            # -3 2 2 2
                            # -3 2 2 5 , in fact, illegal
                            # 
                            if cbk.mon.dihedrallist[j][2] >0: #-+++
                                #
                                if cbk.mon.dihedrallist[j][3] < 0:
                                    print("illegal dihedral")
                                    print("2nd bead >0, 3rd bead <0, illegal")
                                    return -1

                                elif cbk.mon.dihedrallist[j][1]<0 : # this if is redundant.
                                    tmp1bdindex = -cbk.mon.dihedrallist[j][1] - cbk.mon.Natom
                                    tmp1bdnm    = cbk.mon.beadlist[-cbk.mon.dihedrallist[j][1] - 1].name
                                if cbk.mon.dihedrallist[j][2] == cbk.mon.dihedrallist[j][3]:
                                    tmp3bdindex = cbk.mon.dihedrallist[j][3] + cbk.mon.Natom
                                    if cbk.mon.dihedrallist[j][4] == cbk.mon.dihedrallist[j][3]:
                                        tmp4bdindex = cbk.mon.dihedrallist[j][4] + 2* cbk.mon.Natom
                                elif cbk.mon.dihedrallist[j][3] == cbk.mon.dihedrallist[j][4]:
                                    tmp4bdindex = cbk.mon.dihedrallist[j][4] + cbk.mon.Natom
                                    #
                            elif cbk.mon.dihedrallist[j][3] >0: # --++
                                tmp2bdindex = -cbk.mon.dihedrallist[j][2] -    cbk.mon.Natom
                                tmp2bdnm    = cbk.mon.beadlist[ -cbk.mon.anglelist[j][2] -1 ].name
                                if cbk.mon.dihedrallist[j][1] == cbk.mon.dihedrallist[j][2]: #<0
                                    tmp1bdindex = -cbk.mon.dihedrallist[j][1] - 2* cbk.mon.Natom
                                    tmp1bdnm    = cbk.mon.beadlist[ -cbk.mon.anglelist[j][1] -1 ].name
                                else:
                                    tmp1bdindex = -cbk.mon.dihedrallist[j][1] -   cbk.mon.Natom
                                    tmp1bdnm    = cbk.mon.beadlist[ -cbk.mon.anglelist[j][1] -1 ].name
                            elif cbk.mon.dihedrallist[j][4]>0: # ---+
                                tmp3bdindex = -cbk.mon.dihedrallist[j][3] -    cbk.mon.Natom
                                tmp3bdnm    = cbk.mon.beadlist[ -cbk.mon.anglelist[j][3] -1 ].name
                                #
                                if cbk.mon.dihedrallist[j][3] == cbk.mon.dihedrallist[j][2]: #<0
                                    tmp2bdindex = -cbk.mon.dihedrallist[j][2] - 2* cbk.mon.Natom
                                    tmp2bdnm    = cbk.mon.beadlist[ -cbk.mon.anglelist[j][2] -1 ].name
                                    if cbk.mon.dihedrallist[j][1] == cbk.mon.dihedrallist[j][2]: #<0
                                        tmp1bdindex = -cbk.mon.dihedrallist[j][1] - 3* cbk.mon.Natom
                                        tmp1bdnm    = cbk.mon.beadlist[ -cbk.mon.anglelist[j][1] -1 ].name
                                    else:
                                        tmp1bdindex = -cbk.mon.dihedrallist[j][1] -  2* cbk.mon.Natom
                                        tmp1bdnm    = cbk.mon.beadlist[ -cbk.mon.anglelist[j][1] -1 ].name

                                else :
                                    tmp2bdindex = -cbk.mon.dihedrallist[j][2] - cbk.mon.Natom
                                    tmp2bdnm    = cbk.mon.beadlist[ -cbk.mon.anglelist[j][2] -1 ].name
                                    tmp1bdindex = -cbk.mon.dihedrallist[j][1] - cbk.mon.Natom
                                    tmp1bdnm    = cbk.mon.beadlist[ -cbk.mon.anglelist[j][1] -1 ].name
                            # dihedral has more restrict format than angle.
                            # searching for the next monomer is not allowed.
                            # 
                            # In angle def, it can search in the next monomer,
                            # in some special cases.
                            # e.g. [1, -2,2,2] for lmon = 0.0,
                            # but the recommended def should be [1,    -2, -2, 2]
                            # bond need to be revised.
                        # write:
                        if r == 0 and any(bn <= 0 for bn in cbk.mon.dihedrallist[j]):
                            continue
                        elif r ==1:
                            if cbk.mon.dihedrallist[j][2] <0 and cbk.mon.lmon == 0.0:
                                continue

                        elif r ==2:
                            if cbk.mon.anglelist[j][3]<0 and ckb.mon.lmon == 0.0 :
                                continue

                        else: pass
                        self.f.write(\
                        '%8d'  %cbk.mon.dihedrallist[j][0] \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp1bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp2bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp3bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp4bdindex) \
                        + '\n')
                preNatom += cbk.Natom
            self.f.close()
