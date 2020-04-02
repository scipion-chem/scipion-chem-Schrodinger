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
import sys

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import AtomicStructHandler
from schrodinger import Plugin
from schrodinger.objects import SchrodingerMaestroFile

class ProtSchrodingerPrepWizard(EMProtocol):
    """Calls the preparation wizard"""
    _label = 'preparation wizard'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct",
                       label='Atomic Structure:', allowsNull=False)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('preparationStep',self.inputStructure.get().getFileName())
        self._insertFunctionStep('createOutput')

    def preparationStep(self, structFileName):
        prog=Plugin.getHome('utilities/prepwizard')

        fnIn = self._getExtraPath("atomStructIn.pdb")
        aStruct1 = AtomicStructHandler(structFileName)
        aStruct1.write(fnIn)

        args='-WAIT extra/atomStructIn.pdb atomStructOut.mae'
        self.runJob(prog,args,cwd=self._getPath())

    def createOutput(self):
        fnMae = self._getPath('atomStructOut.mae')
        if os.path.exists(fnMae):
            maeFile=SchrodingerMaestroFile()
            maeFile.setFileName(fnMae)

            self._defineOutputs(outputMae=maeFile)
            self._defineSourceRelation(self.inputStructure, maeFile)
