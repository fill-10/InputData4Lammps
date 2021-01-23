# class_polymer.py
#
#
from class_block import block

class polymer:
    def __init__(self, blocklist):
        # create an empty list of blocks
        self.block_list = []
        self.TotNatom = 0
        # add blocks to the list
        #i = 0
        for blk in blocklist:
            #if i == 0: # first
            #    blk.prebond = None
            #if i == len(blocklist)-1: # last
            #    blk.sucbond = None
            self.block_list += [blk]
            # calculate the total number of beads in the polymer
            self.TotNatom += blk.Natom

    def writexyz(self, filename='chain1.xyz'):
        f= open(filename, 'w')
        f.write(str(self.TotNatom)+'\n\n')
        pre_block_Vec = [0, 0, 0]
        # 
        iblk = 0
        for blk in self.block_list:
            pre_block_Vec[0] += blk.prelinkvector[0]
            pre_block_Vec[1] += blk.prelinkvector[1]
            pre_block_Vec[2] += blk.prelinkvector[2]
            for j in range(0, blk.monrep):
                for bd in blk.mon.beadlist:
                    f.write( '%s' %bd.name \
                    + 3*' ' \
                    + '%6.2f' %( pre_block_Vec[0] + bd.coor[0] \
                    + j * ( blk.mon.deltaV[0] + blk.monlink[0] )  )\
                    + 3*' ' \
                    + '%6.2f' %( pre_block_Vec[1] + bd.coor[1] \
                    + j * ( blk.mon.deltaV[1] + blk.monlink[1] )  )\
                    + 3*' ' \
                    + '%6.2f' %( pre_block_Vec[2] + bd.coor[2] \
                    + j * ( blk.mon.deltaV[2] + blk.monlink[2] )  )\
                    + 3*' ' \
                    + '\n' \
                    )
            iblk += 1
            pre_block_Vec[0] += blk.monrep    * blk.mon.deltaV[0] \
                                + (blk.monrep-1) * blk.monlink[0]
            pre_block_Vec[1] += blk.monrep    * blk.mon.deltaV[1] \
                                + (blk.monrep-1) * blk.monlink[1]
            pre_block_Vec[2] += blk.monrep    * blk.mon.deltaV[2] \
                                + (blk.monrep-1) * blk.monlink[2]

        f.close()
	#
	#
    def writebond(self, filename):
        f = open(filename, 'w+')
        preNatom= 0
        icbk = 0
        L_prebond = [ ]

        for cbk in self.block_list:
            L_prebond = cbk.prebond # pre bond list
            for pbnd in L_prebond:
                f.write( \
                '%8d' %pbnd[0] \
                + '%12d' %(  preNatom + pbnd[1] +1 ) \
                + '%12d' %(  preNatom + pbnd[2] +1 ) \
                + '%8s' %pbk.mon.beadlist[ pbnd[1]%pbk.mon.Natom ].name \
                + '%8s' %cbk.mon.beadlist[ pbnd[2]%cbk.mon.Natom ].name \
                + '\n')

            for r in range(0, cbk.monrep):
                # bonds inside one monomer (normal)
                for bnd in cbk.mon.bondlist:
                    # index starts from 0 in python
                    bd1name = cbk.mon.beadlist[ bnd[1]%cbk.mon.Natom ].name
                    bd2name = cbk.mon.beadlist[ bnd[2]%cbk.mon.Natom ].name
                    if bnd[1] < - r * cbk.mon.Natom :
                        continue #skip the link bond if the first mon
                    f.write(\
                    '%8d' %bnd[0] \
                    + '%12d'   %(preNatom + r * cbk.mon.Natom + bnd[1]+1) \
                    + '%12d'   %(preNatom + r * cbk.mon.Natom + bnd[2]+1) \
                    + '%8s'    %bd1name \
                    + '%8s'    %bd2name \
                    + '\n')
            #
            #
            preNatom += cbk.Natom
            pbk = self.block_list[icbk] # pre block 
            icbk += 1

        f.close()
	#
	#
	#
    def writeangle(self, filename):
        f=open(filename, 'w')
        preNatom= 0
        icbk = 0
        L_preang = []
        for cbk in self.block_list:
            L_preang = cbk.headangle_list
            for preang in L_preang: # pre angles
                bd1name = pbk.mon.beadlist[ preang[1]%pbk.mon.Natom ].name
                # this code cannot deal with the angle spanning over 3 blocks
                if preang[2]>=0:
                    bd2name = cbk.mon.beadlist[ preang[2]%cbk.mon.Natom ].name
                else:
                    bd2name = pbk.mon.beadlist[ preang[2]%pbk.mon.Natom ].name
                bd3name = cbk.mon.beadlist[preang[3]%cbk.mon.Natom].name
                f.write( \
                '%8d'    %preang[0] \
                + '%12d' %(  preNatom + preang[1] +1 ) \
                + '%12d' %(  preNatom + preang[2] +1 ) \
                + '%12d' %(  preNatom + preang[3] +1 ) \
                + '%8s'  %bd1name \
                + '%8s'  %bd2name \
                + '%8s'  %bd3name \
                + '\n' )
            for r in range(0, cbk.monrep): # inside block
                for angle in cbk.mon.anglelist:
                    if angle[1] < -r * cbk.mon.Natom \
                    or angle[2] < -r * cbk.mon.Natom :
                        continue
                    f.write( \
                    '%8d'  %angle[0] \
                    + '%12d' %( preNatom + r * cbk.mon.Natom + angle[1] +1) \
                    + '%12d' %( preNatom + r * cbk.mon.Natom + angle[2] +1) \
                    + '%12d' %( preNatom + r * cbk.mon.Natom + angle[3] +1) \
                    + '%8s'  %cbk.mon.beadlist[angle[1]%cbk.mon.Natom ].name\
                    + '%8s'  %cbk.mon.beadlist[angle[2]%cbk.mon.Natom ].name\
                    + '%8s'  %cbk.mon.beadlist[angle[3]%cbk.mon.Natom ].name\
                    + '\n')

            ##
            preNatom += cbk.Natom
            pbk = self.block_list[icbk]
            icbk += 1

        f.close()
       
    def writedihedral(self, filename):
        f= open(filename, 'w')
        preNatom= 0 
        icbk = 0
        L_predih = []
        for cbk in self.block_list:
            L_predih = cbk.headdih_list
            for predih in L_predih:
                bd1name = pbk.mon.beadlist[ predih[1]%pbk.mon.Natom ].name
                if predih[2]>=0:
                    bd2name = cbk.mon.beadlist[ predih[2]%cbk.mon.Natom ].name
                else:
                    bd2name = pbk.mon.beadlist[ predih[2]%pbk.mon.Natom ].name
                if predih[3]>=0:
                    bd3name = cbk.mon.beadlist[ predih[3]%cbk.mon.Natom ].name
                else:
                    bd3name = pbk.mon.beadlist[ predih[3]%pbk.mon.Natom ].name
                bd4name = cbk.mon.beadlist[predih[4]%cbk.mon.Natom].name
                f.write( \
                '%8d'    %predih[0] \
                + '%12d' %(  preNatom + predih[1] +1 ) \
                + '%12d' %(  preNatom + predih[2] +1 ) \
                + '%12d' %(  preNatom + predih[3] +1 ) \
                + '%12d' %(  preNatom + predih[4] +1 ) \
                + '%8s'  %bd1name \
                + '%8s'  %bd2name \
                + '%8s'  %bd3name \
                + '%8s'  %bd4name \
                + '\n' )
            for r in range(0, cbk.monrep): # inside block
                for dih in cbk.mon.dihedrallist:
                    if dih[1] < -r * cbk.mon.Natom \
                    or dih[2] < -r * cbk.mon.Natom \
                    or dih[3] < -r * cbk.mon.Natom :
                        continue
                    f.write( \
                    '%8d'  %dih[0] \
                    + '%12d' %( preNatom + r * cbk.mon.Natom + dih[1] +1) \
                    + '%12d' %( preNatom + r * cbk.mon.Natom + dih[2] +1) \
                    + '%12d' %( preNatom + r * cbk.mon.Natom + dih[3] +1) \
                    + '%12d' %( preNatom + r * cbk.mon.Natom + dih[4] +1) \
                    + '%8s'  %cbk.mon.beadlist[dih[1]%cbk.mon.Natom ].name\
                    + '%8s'  %cbk.mon.beadlist[dih[2]%cbk.mon.Natom ].name\
                    + '%8s'  %cbk.mon.beadlist[dih[3]%cbk.mon.Natom ].name\
                    + '%8s'  %cbk.mon.beadlist[dih[4]%cbk.mon.Natom ].name\
                    + '\n')
            ##
            preNatom += cbk.Natom
            pbk = self.block_list[icbk]
            icbk += 1

        f.close()
                
    def writeimproper(self, filename):
        f= open(filename, 'w')
        preNatom= 0 
        icbk = 0
        L_preimp = []
        for cbk in self.block_list:
            L_preimp = cbk.headimproper_list
            for preimp in L_preimp:
                bd1name = pbk.mon.beadlist[ preimp[1]%pbk.mon.Natom ].name
                if preimp[2]>=0:
                    bd2name = cbk.mon.beadlist[ preimp[2]%cbk.mon.Natom ].name
                else:
                    bd2name = pbk.mon.beadlist[ preimp[2]%pbk.mon.Natom ].name
                if preimp[3]>=0:
                    bd3name = cbk.mon.beadlist[ preimp[3]%cbk.mon.Natom ].name
                else:
                    bd3name = pbk.mon.beadlist[ preimp[3]%pbk.mon.Natom ].name
                bd4name = cbk.mon.beadlist[preimp[4]%cbk.mon.Natom].name
                f.write( \
                '%8d'    %preimp[0] \
                + '%12d' %(  preNatom + preimp[1] +1 ) \
                + '%12d' %(  preNatom + preimp[2] +1 ) \
                + '%12d' %(  preNatom + preimp[3] +1 ) \
                + '%12d' %(  preNatom + preimp[4] +1 ) \
                + '%8s'  %bd1name \
                + '%8s'  %bd2name \
                + '%8s'  %bd3name \
                + '%8s'  %bd4name \
                + '\n' )
            for r in range(0, cbk.monrep): # inside block
                for imp in cbk.mon.improperlist:
                    if imp[1] < -r * cbk.mon.Natom \
                    or imp[2] < -r * cbk.mon.Natom \
                    or imp[3] < -r * cbk.mon.Natom :
                        continue
                    f.write( \
                    '%8d'  %imp[0] \
                    + '%12d' %( preNatom + r * cbk.mon.Natom + imp[1] +1) \
                    + '%12d' %( preNatom + r * cbk.mon.Natom + imp[2] +1) \
                    + '%12d' %( preNatom + r * cbk.mon.Natom + imp[3] +1) \
                    + '%12d' %( preNatom + r * cbk.mon.Natom + imp[4] +1) \
                    + '%8s'  %cbk.mon.beadlist[imp[1]%cbk.mon.Natom ].name\
                    + '%8s'  %cbk.mon.beadlist[imp[2]%cbk.mon.Natom ].name\
                    + '%8s'  %cbk.mon.beadlist[imp[3]%cbk.mon.Natom ].name\
                    + '%8s'  %cbk.mon.beadlist[imp[4]%cbk.mon.Natom ].name\
                    + '\n')
            ##
            preNatom += cbk.Natom
            pbk = self.block_list[icbk]
            icbk += 1

        f.close()

