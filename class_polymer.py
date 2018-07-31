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
        def writebond(self,filename):
                self.f = open(filename, 'w+')
                preNatom= 0 #set up a counter to count the number of atoms in all previously covered blocks.
                #
                for i in range(0, len(self.block_list)):
                        cbk = self.block_list[i] # current block
			#
                        # bond linking the previous blocks, starting from the second block:
                        if i>0:
                                pbk = self.block_list[i-1] # read the previous block
                                if cbk.prebond:
                                    # check:
                                    if cbk.prebond[1] >0:
                                        print("illegal block prebond:")
                                        print(i)
                                        return -1
                                    #
                                    tmp1bdindex = cbk.prebond[1]
                                    tmp2bdindex = cbk.prebond[2]
                                    tmp1bdname  = pbk.mon.beadlist[tmp1bdindex].name
                                    tmp2bdname  = cbk.mon.beadlist[tmp2bdindex-1].name

                                    self.f.write(\
                                    '%8d' %cbk.prebond[0] \
                                    + '%12d' %(  preNatom + tmp1bdindex +1 ) \
                                    + '%12d' %(  preNatom + tmp2bdindex    ) \
                                    + '%8s' %tmp1bdname \
                                    + '%8s' %tmp2bdname \
                                    + '\n')


			for r in range(0, cbk.monrepeat):
				# bonds inside one monomer (normal)
				for j in range(0, len(cbk.mon.bondlist)):
                                    # regular bond inside monomer
                                    tmp1bdindex = cbk.mon.bondlist[j][1]
                                    tmp2bdindex = cbk.mon.bondlist[j][2]
                                    #
                                    if tmp1bdindex > 0:
                                        tmp1bdname  = cbk.mon.beadlist[tmp1bdindex-1].name
                                    else: # when bead1 is negative
                                        tmp1bdname  = cbk.mon.beadlist[tmp1bdindex].name
                                    #
                                    tmp2bdname  = cbk.mon.beadlist[tmp2bdindex-1].name
                                    #
                                    if cbk.mon.bondlist[j][1] <0 and r ==0:
                                        continue # skip the linking bond for the first monomer r = 0

                                    self.f.write(\
                                    '%8d' %cbk.mon.bondlist[j][0] \
                                    + '%12d'   %(preNatom + r * cbk.mon.Natom + tmp1bdindex) \
                                    + '%12d'   %(preNatom + r * cbk.mon.Natom + tmp2bdindex) \
                                    + '%8s'    %tmp1bdname \
                                    + '%8s' %tmp2bdname\
                                    + '\n')

			preNatom += cbk.Natom

		self.f.close()
	#
	#
	#
	def writeangle(self,filename):
            self.f=open(filename, 'w+')
            preNatom= 0 #set up a counter to count the number of atoms in all previous blocks.
            for i in range(0, len(self.block_list)):
                cbk = self.block_list[i] # current block
		#
		if i>0:    # Write the angle linking the previous block.
                    pbk = self.block_list[i-1] # read the previous block
                    if cbk.headangle_list:
                        for j in range(0, len(cbk.headangle_list) ):
                            #
                            tmp1bdindex = cbk.headangle_list[j][1]
                            tmp2bdindex = cbk.headangle_list[j][2]
                            tmp3bdindex = cbk.headangle_list[j][3]
                            tmp3bdname =  cbk.mon.beadlist[tmp3bdindex-1].name
                            # if the apex is in cbk
                            if cbk.headangle_list[j][2]>0 :
                                #
                                tmp2bdname = cbk.mon.beadlist[tmp2bdindex - 1].name
                                tmp1bdname = pbk.mon.beadlist[tmp1bdindex].name
                                tmp1bdindex += 1
                                #
                                # special case: -4 1 1  and -1 2 2
                                # very special
                                if cbk.headangle_list[j][2] == cbk.headangle_list[j][3]:
                                    tmp3bdindex += cbk.mon.Natom
                            #
                            # if apex bead in pbk:
                            if cbk.headangle_list[j][2]< 0 :
                                tmp2bdname = pbk.mon.beadlist[tmp2bdindex].name
                                tmp2bdindex += 1
                                tmp1bdname = pbk.mon.beadlist[tmp1bdindex%pbk.mon.Natom].name
                                tmp1bdindex += 1
                            self.f.write(\
                            '%8d'    %cbk.headangle_list[j][0] \
                            + '%12d'  %(  preNatom + tmp1bdindex  ) \
                            + '%12d' %(  preNatom + tmp2bdindex  ) \
                            + '%12d' %(  preNatom + tmp3bdindex  ) \
                            + '%8s'  %tmp1bdname \
                            + '%8s'  %tmp2bdname \
                            + '%8s'  %tmp3bdname \
                            + '\n' )
                        #
                        try: del j
                        except: pass


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
                #if len(cbk.mon.anglelist) == 0:
                #    preNatom += cbk.Natom
                #    continue
                #
                #
                ###

                for r in range(0, cbk.monrepeat):
                # From the second monomer, need to include the angles linking the current and previous monomers.
                # I.e., r from 1 to self.repeat-1 (incld)
                    for j in range(0, len(cbk.mon.anglelist)):
                        try:
                            tmp1bdindex = cbk.mon.anglelist[j][1] 
                            tmp1bdname = cbk.mon.beadlist[tmp1bdindex-1].name
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp2bdindex = cbk.mon.anglelist[j][2] 
                            tmp2bdname = cbk.mon.beadlist[tmp2bdindex-1].name
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp3bdindex = cbk.mon.anglelist[j][3] 
                            tmp3bdname = cbk.mon.beadlist[tmp3bdindex-1].name
                        except:
                            pass
                        finally:
                            pass
                        #
                        # if the bond is a monomer linking bond,  update:
                        #
                        if cbk.mon.anglelist[j][2] >0 and cbk.mon.anglelist[j][1] <0:
                            tmp1bdname  =  cbk.mon.beadlist[ tmp1bdindex ].name
                            tmp1bdindex += 1
                        if cbk.mon.anglelist[j][2] <0 and cbk.mon.anglelist[j][1] <0:
                            tmp1bdname  =  cbk.mon.beadlist[tmp1bdindex%cbk.mon.Natom].name
                            tmp1bdindex += 1
                            tmp2bdname  =  cbk.mon.beadlist[tmp2bdindex%cbk.mon.Natom ].name
                            tmp2bdindex += 1

                        #write:
                        if r==0 and any(bn <= 0 for bn in cbk.mon.anglelist[j]):
                            # if it is the monomer linking angle
                            # and this is the first monomer
                            # go to the next angle, i.e. j++
                            continue
                        
                        elif r ==1:
                            # second monomer
                            if cbk.mon.anglelist[j][2]<0 and cbk.mon.anglelist[j][1]<0 \
                            and cbk.mon.anglelist[j][1]/cbk.mon.Natom < -1:
                                continue
                        #print(tmp1bdindex, tmp2bdindex, tmp3bdindex, tmp1bdnm, tmp2bdnm,tmp3bdnm)
                        self.f.write(\
                        '%8d'  %cbk.mon.anglelist[j][0] \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp1bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp2bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp3bdindex) \
                        + '%8s'  %tmp1bdname \
                        + '%8s'  %tmp2bdname \
                        + '%8s'  %tmp3bdname \
                        + '\n')
                #                
                preNatom += cbk.Natom
		
            self.f.close()
        
        def writedihedral(self, filename):
            self.f= open(filename, 'w+')
            preNatom= 0 #set up a counter to count the number of atoms in all previously covered blocks.
            for i in range(0, len(self.block_list)):
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
                                tmp2bdindex += 1
                                tmp3bdindex += 1
                                tmp1bdindex += 1
                            
                            elif cbk.headdih_list[m][2]<0:
                                tmp1bdindex += 1
                                tmp2bdindex += 1
                                if cbk.headdih_list[m][3] == cbk.headdih_list[m][4]:
                                    tmp4bdindex += cbk.mon.Natom

                            elif cbk.headdih_list[m][1]<0:
                                tmp1bdindex += 1
                                if ckb.headdih_list[m][2] == ckb.headdih_list[m][3]:
                                    tmp3bdidnex += cbk.mon.Natom
                                    if cbk.headdih_list[m][3] == cbk.headdih_list[m][4]:
                                        tmp4bdindex += cbk.mon.Natom

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
                            #tmp1bdnm = cbk.mon.beadlist[tmp1bdindex-1].name
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp2bdindex = cbk.mon.dihedrallist[j][2] 
                            #tmp2bdnm = cbk.mon.beadlist[tmp2bdindex-1].name
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp3bdindex = cbk.mon.dihedrallist[j][3] 
                            #tmp3bdnm = cbk.mon.beadlist[tmp3bdindex-1].name 
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp4bdindex = cbk.mon.dihedrallist[j][4] 
                            #tmp4bdnm = cbk.mon.beadlist[tmp3bdindex-1].name 
                        except:
                            pass
                        finally:
                            pass

                        #update temp bd index for negative beads
                        if any(bn <= 0 for bn in cbk.mon.dihedrallist[j]):
                            #
                            # 
                            if cbk.mon.dihedrallist[j][2] >0: #-+++

                                if cbk.mon.dihedrallist[j][1]<0 : # this if is redundant.
                                    tmp1bdnm += 1
                                    #
                            elif cbk.mon.dihedrallist[j][3] >0: # --++
                                tmp2bdindex += 1
                                tmp1bdindex += 1
                            elif cbk.mon.dihedrallist[j][4]>0: # ---+\
                                tmp3bdindex += 1
                                tmp2bdindex += 1
                                tmp1bdindex += 1
                        #
                        #
                        # write:
                        if r == 0 and any(bn <= 0 for bn in cbk.mon.dihedrallist[j]):
                            continue
                        if r ==1 and any(bn/cbk.mon.Natom < -1 for bn in cbk.mon.dihedrallist[j]): #--++
                            continue
                        if r ==2 and any(bn/cbk.mon.Natom < -2 for bn in cbk.mon.dihedrallist[j]): # ---+
                            continue

                        self.f.write(\
                        '%8d'  %cbk.mon.dihedrallist[j][0] \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp1bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp2bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp3bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp4bdindex) \
                        + '\n')
                preNatom += cbk.Natom
            self.f.close()
       
        def writeimproper(self, filename):
            self.f= open(filename, 'w+')
            preNatom= 0 #set up a counter to count the number of atoms in all previously covered blocks.
            for i in range(0, len(self.block_list)):
                ###--------------------------###
                ### block linking improper   ###
                ###--------------------------###
                cbk = self.block_list[i]
                if i > 0:
                    pbk = self.block_list[i-1]
                    if cbk.headimproper_list:
                        for m in range(0, len(cbk.headimproper_list) ):
                            tmp1bdindex = cbk.headimproper_list[m][1]
                            tmp2bdindex = cbk.headimproper_list[m][2]
                            tmp3bdindex = cbk.headimproper_list[m][3]
                            tmp4bdindex = cbk.headimproper_list[m][4]
                            # translate
                            if cbk.headimproper_list[m][3]<0: #---+
                                tmp2bdindex += 1
                                tmp3bdindex += 1
                                tmp1bdindex += 1
                            
                            elif cbk.headimproper_list[m][2]<0: #--++
                                tmp1bdindex += 1
                                tmp2bdindex += 1
                                if cbk.headimproper_list[m][3] == cbk.headimproper_list[m][4]:
                                    tmp4bdindex += cbk.mon.Natom

                            elif cbk.headimproper_list[m][1]<0: #-+++
                                tmp1bdindex += 1
                                if ckb.headimproper_list[m][2] == ckb.headimproper_list[m][3]:
                                    tmp3bdidnex += cbk.mon.Natom
                                    if cbk.headimproper_list[m][3] == cbk.headimproper_list[m][4]:
                                        tmp4bdindex += cbk.mon.Natom

                            # write
                            self.f.write(\
                            '%8d'  %cbk.headimproper_list[m][0] \
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
                ### improper in each block   ###
                ###--------------------------###
                ###
                # This IF is a must.
                # There is the case:
                # only block.headimproper_list has improper,
                # but no improper inside each block,
                # then the headimproper will pollute the following lines
                # and cause error.
                if len(cbk.mon.improperlist) == 0:
                    preNatom += cbk.Natom
                    continue
                #
                #
                ###
                for r in range(0, cbk.monrepeat):
                    for j in range(0 , len(cbk.mon.improperlist)):
                        try:
                            tmp1bdindex = cbk.mon.improperlist[j][1] 
                            #tmp1bdnm = cbk.mon.beadlist[tmp1bdindex-1].name
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp2bdindex = cbk.mon.improperlist[j][2] 
                            #tmp2bdnm = cbk.mon.beadlist[tmp2bdindex-1].name
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp3bdindex = cbk.mon.improperlist[j][3] 
                            #tmp3bdnm = cbk.mon.beadlist[tmp3bdindex-1].name 
                        except:
                            pass
                        finally:
                            pass
                        try:
                            tmp4bdindex = cbk.mon.improperlist[j][4] 
                            #tmp4bdnm = cbk.mon.beadlist[tmp3bdindex-1].name 
                        except:
                            pass
                        finally:
                            pass

                        #update temp bd index for negative beads
                        if any(bn <= 0 for bn in cbk.mon.improperlist[j]):
                            #
                            # 
                            if cbk.mon.improperlist[j][2] >0: #-+++

                                if cbk.mon.improperlist[j][1]<0 : # this if is redundant.
                                    tmp1bdnm += 1
                                    #
                            elif cbk.mon.improperlist[j][3] >0: # --++
                                tmp2bdindex += 1
                                tmp1bdindex += 1
                            elif cbk.mon.improperlist[j][4]>0: # ---+\
                                tmp3bdindex += 1
                                tmp2bdindex += 1
                                tmp1bdindex += 1
                        #
                        #
                        # write:
                        if r == 0 and any(bn <= 0 for bn in cbk.mon.improperlist[j]):
                            continue
                        if r ==1 and any(bn/cbk.mon.Natom < -1 for bn in cbk.mon.improperlist[j]): #--++
                            continue
                        if r ==2 and any(bn/cbk.mon.Natom < -2 for bn in cbk.mon.improperlist[j]): # ---+
                            continue

                        self.f.write(\
                        '%8d'  %cbk.mon.improperlist[j][0] \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp1bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp2bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp3bdindex) \
                        + '%12d' %( preNatom + r * cbk.mon.Natom + tmp4bdindex) \
                        + '\n')
                preNatom += cbk.Natom
            self.f.close()


