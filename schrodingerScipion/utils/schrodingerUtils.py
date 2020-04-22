# -*- coding: utf-8 -*-
#  **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import sys
import schrodinger.structure as structure

def getCentroid(st):
    """
    Return the (x,y,z) coordinates of the centroid of the structure 'st'.
    """
    natoms = len(st.atom)
    if natoms == 0:
        raise Exception("Error: No atoms in structure. Can't calculate centroid.")
    xsum = 0.0
    ysum = 0.0
    zsum = 0.0
    for at in st.atom:
        xsum += at.x
        ysum += at.y
        zsum += at.z
    return (xsum / natoms, ysum / natoms, zsum / natoms)

def readStructure(fnIn):
    if "@" in fnIn:
        number, fnIn = fnIn.split('@')
        number = int(number)
    else:
        number = 0
    nst = structure.count_structures(fnIn)
    if number>nst:
        return None
    else:
        sr = structure.StructureReader(fnIn)
        i = 0
        for st in sr:
            if i==number:
                return st
            i+=1

if __name__ == "__main__":
    if len(sys.argv)==1:
        print("Usage: python3 schrodingerUtils.py [options]")
        print("   centroid <structureFile> <centroid.txt>")
    elif sys.argv[1]=="centroid":
        fnIn = sys.argv[2]
        fnOut = sys.argv[3]
        st = readStructure(fnIn)
        x,y,z = getCentroid(st)
        fh=open(fnOut,'w')
        fh.write("%f %f %f\n"%(x,y,z))
        fh.close()
