# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from pyworkflow.utils.path import moveFile

def putMol2Title(fn, title=""):
    i=0
    fhIn = open(fn)
    fhOut = open(fn+".aux",'w')
    for line in fhIn.readlines():
        if i!=1:
            fhOut.write(line)
        else:
            if title!="":
                fhOut.write(title+"\n")
            else:
                fhOut.write(os.path.splitext(os.path.split(fn)[1])[0]+"\n")
        i+=1
    fhIn.close()
    fhOut.close()
    moveFile(fn+".aux",fn)