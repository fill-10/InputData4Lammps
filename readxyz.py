# this file is to convert xyz file into bead instance.
# TODO: under development
import os
from class_bead import bead
def readxyz(xyzfilename):
    xyz_f = open(xyzfilename, 'r')
    xyz_f.seek(0)
    rawdata = xyz_f.readlines()
    nmol = int(rawdata[0].split()[0])
    print nmol
    tmpbead = bead()
