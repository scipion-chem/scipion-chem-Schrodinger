# **************************************************************************
# *
# * Authors: Martín Salinas Antón (martin.salinas@cnb.csic.es)
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

from pwem.protocols import EMProtocol
from schrodingerScipion import Plugin

class ProtSchrodingerMMGBSA(EMProtocol):
    """ TODO: find out & place here"""
    _label = 'mm-gbsa'
    _program = ""

    def _defineParams(self, form):
        pass

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        """ This function inserts all steps functions that will be run when running the protocol """
        self._insertFunctionStep('runMMGBSAStep')
    
    def runMMGBSAStep(self):
        """ This function runs the schrodinger binary file with the given params """
        # Getting binary file
        schrodinger = self.getMMGBSABinaryFile()
        print("BINARY FILE:", schrodinger)
    
    # --------------------------- Utils functions --------------------
    def getMMGBSABinaryFile(self):
        """ This function returns the location for the Schrodinger MM-GBSA binary file """
        return os.path.join(Plugin.getVar('SCHRODINGER_HOME'), 'prime_mmgbsa')