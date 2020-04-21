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
import glob
import os

from pyworkflow.protocol.params import PointerParam
from pyworkflow.utils.path import createLink

from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct
from pwem.convert.atom_struct import AtomicStructHandler

from schrodinger import Plugin as schrodinger_plugin
from schrodinger.objects import SchrodingerGrid

class ProtSchrodingerGridManual(EMProtocol):
    """Calls glide GUI to prepare a grid"""
    _label = 'grid definition manual (glide)'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct, SchrodingerAtomStruct",
                       label='Atomic Structure:', allowsNull=False)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('preparationStep')
        self._insertFunctionStep('createOutput')

    def preparationStep(self):
        if isinstance(self.inputStructure.get(),AtomStruct):
            fnIn = self._getExtraPath("atomStructIn.pdb")
            aStruct1 = AtomicStructHandler(self.inputStructure.get().getFileName())
            aStruct1.write(fnIn)
            fnIn='extra/atomStructIn.pdb'
            self.runJob(schrodinger_plugin.getHome('maestro'), "-b %s"%fnIn, cwd=self._getPath())
        else:
            fnIn = self._getExtraPath("atomStructIn")+self.inputStructure.get().getExtension()
            createLink(self.inputStructure.get().getFileName(),fnIn)
            fnIn=os.path.join('extra',os.path.split(fnIn)[1])
            self.runJob(schrodinger_plugin.getHome('maestro'), "-m %s"%fnIn, cwd=self._getPath())

    def createOutput(self):
        for fnDir in glob.glob(self._getPath('glide-*')):
            fnBase = os.path.split(fnDir)[1]
            fnGrid = os.path.join(fnDir,'%s.zip'%fnBase)
            if os.path.exists(fnGrid):
                gridFile=SchrodingerGrid(filename=fnGrid)
                gridFile.setStructure(self.inputStructure)
                self._defineOutputs(outputGrid=gridFile)
                self._defineSourceRelation(self.inputStructure, gridFile)
