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
from pyworkflow.object import String

from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct
from pwem.convert.atom_struct import AtomicStructHandler

from schrodingerScipion import Plugin as schrodinger_plugin
from schrodingerScipion.objects import SchrodingerGrid, SetOfSchrodingerGrids

class ProtSchrodingerGridManual(EMProtocol):
    """Calls glide GUI to prepare a grid.
    Please do not change the default JobNames and close Maestro only after all the jobs are finished"""
    _label = 'grid definition manual (glide)'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct", label='Atomic Structure:',
                      help='Input structure to generate the schrodinger grids on.'
                           'A tutorial for grid generation can be found at '
                           'https://www.youtube.com/watch?v=_AUKLGtrBR8')

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
        outGrids = SetOfSchrodingerGrids(filename=self._getPath('SchGrids.sqlite'))
        fnStruct = glob.glob(self._getExtraPath("atomStructIn*"))[0]

        for fnDir in glob.glob(self._getPath('glide-*')):
            fnBase = os.path.split(fnDir)[1]
            gridFile, gridArgs = self.getGridArgs(fnBase)

            fnGrid = os.path.join(fnDir, gridFile)
            gridId = int(fnDir.split('glide-grid_')[1])
            if os.path.exists(fnGrid):
                gridObj = SchrodingerGrid(filename=fnGrid, **gridArgs)
                gridObj._proteinFile = String(fnStruct)
                gridObj.setObjId(gridId)

                oriFile = self.inputStructure.get().getFileName()
                if '.mae' in oriFile:
                    gridObj.structureFile = String(oriFile)

                outGrids.append(gridObj)

        outGrids.buildBBoxesPML()
        self._defineOutputs(outputGrids=outGrids)
        self._defineSourceRelation(self.inputStructure, outGrids)


    def getGridArgs(self, fnBase):
        with open(self._getPath('{}/{}.in'.format(fnBase, fnBase))) as f:
            for line in f:
                if line.startswith('FORCEFIELD'):
                    ff = line.split()[1].strip()

                elif line.startswith('GRID_CENTER'):
                    cMass = line.replace(',', '').split()[1:]

                elif line.startswith('GRIDFILE'):
                    gridFile = line.split()[1].strip()

                elif line.startswith('INNERBOX'):
                    iBox = line.replace(',', '').split()[1:]

                elif line.startswith('OUTERBOX'):
                    oBox = line.replace(',', '').split()[1:]

        return gridFile, {'centerX': cMass[0], 'centerY': cMass[1], 'centerZ': cMass[2],
                          'innerX': iBox[0], 'innerY': iBox[1], 'innerZ': iBox[2],
                          'outerX': oBox[0], 'outerY': oBox[1], 'outerZ': oBox[2],
                          'forceField': ff}
