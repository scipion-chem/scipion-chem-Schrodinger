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

from pyworkflow.protocol.params import PointerParam
from pyworkflow.utils.path import createLink

from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct
from pwem.convert.atom_struct import AtomicStructHandler

from schrodingerScipion import Plugin as schrodinger_plugin
from schrodingerScipion.objects import SchrodingerAtomStruct

class ProtSchrodingerPrepWizardManual(EMProtocol):
    """Calls the preparation wizard GUI.
    Please do not change the default JobNames and close Maestro only after all the jobs are finished"""
    _label = 'target preparation manual (prepwizard)'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct", label='Atomic Structure to prepare: ',
                      help='Input structure to prepare. A tutorial can be found at '
                           'https://www.youtube.com/watch?v=YRFROyN88Fw&ab_channel=Schr%C3%B6dingerTV')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('preparationStep')
        self._insertFunctionStep('createOutput')

    def preparationStep(self):

        oriFile = self.inputStructure.get().getFileName()
        _, inExt = os.path.splitext(oriFile)

        if inExt in ['.pdb', '.mae', '.maegz']:
            fnIn = os.path.abspath(self._getExtraPath("atomStructIn{}".format(inExt)))
            createLink(oriFile, fnIn)

        else:
            fnIn = os.path.abspath(self._getExtraPath("atomStructIn.pdb"))
            aStruct1 = AtomicStructHandler(oriFile)
            aStruct1.write(fnIn)

        self.runJob(schrodinger_plugin.getHome('maestro'), " %s" % fnIn, cwd=self._getPath())

    def createOutput(self):
        files = [self._getPath('prepwizard_workdir/%s'%file)
                 for file in os.listdir(self._getPath("prepwizard_workdir")) if (file.lower().endswith('out.mae'))]
        if len(files) > 0:
            files.sort(key=os.path.getmtime)
            filesSorted = sorted(files, key=os.path.getmtime)

            maeFile = SchrodingerAtomStruct()
            maeFile.setFileName(filesSorted[-1])

            self._defineOutputs(outputStructure=maeFile)
            self._defineSourceRelation(self.inputStructure, maeFile)

    def _citations(self):
        return ['Sastry2013']
