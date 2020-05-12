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

from pyworkflow.protocol.params import PointerParam, EnumParam, FloatParam, BooleanParam, IntParam
from pyworkflow.utils.path import copyFile, moveFile, cleanPath
from pwem.protocols import EMProtocol
from schrodingerScipion import Plugin
from bioinformatics.objects import SetOfDatabaseID, SetOfSmallMolecules, SmallMolecule

class ProtSchrodingerPrepareForWeb(EMProtocol):
    """Prepare a set of small molecules for the web"""
    _label = 'prepare for web'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSmallMols', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Set of small molecules:', allowsNull=False)

        # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('prepareStep')

    def prepareStep(self):
        def addParam(headerLine, lineToPrint, lineNo, param, value=None):
            if value is None:
                if hasattr(mol,param):
                    value = getattr(mol,param).get()
            if value is not None:
                if lineNo==0:
                    if headerLine!="":
                        headerLine+="; "
                    headerLine+=param
                if lineToPrint!="":
                    lineToPrint+="; "
                lineToPrint += str(value)
            return headerLine, lineToPrint

        def getGridRun(mol):
            fn = mol.poseFile.get()
            return fn.split('Runs/')[1][0:6]

        progStructConvert=Plugin.getHome('utilities/structconvert')
        progMaeSubset=Plugin.getHome('utilities/maesubset')
        grids = {}
        gridNo = 1
        for mol in self.inputSmallMols.get():
            gridRun = getGridRun(mol)
            if not gridRun in grids:
                grids[gridRun] = gridNo
                gridNo += 1

        fh = open(self._getExtraPath("smallMolecules.csv"),'w')
        lineNo = 0
        headerLine = ""
        for mol in self.inputSmallMols.get():
            fnSmall = mol.smallMoleculeFile.get()
            fnMol = os.path.split(fnSmall)[1]
            fnRoot = os.path.splitext(fnMol)[0]
            if "-" in fnRoot:
                fnRoot = fnRoot.split('-')[0]

            lineToPrint = ""
            headerLine, lineToPrint = addParam(headerLine, lineToPrint, lineNo, "smallMolecule", fnRoot)
            if fnRoot.startswith("ZINC"):
                headerLine, lineToPrint = addParam(headerLine, lineToPrint, lineNo, "url",
                                                   "https://zinc15.docking.org/substances/%s"%fnRoot)

            headerLine, lineToPrint = addParam(headerLine, lineToPrint, lineNo, "bindingSiteScore")
            headerLine, lineToPrint = addParam(headerLine, lineToPrint, lineNo, "dockingScore")
            headerLine, lineToPrint = addParam(headerLine, lineToPrint, lineNo, "ligandEfficiency")

            headerLine, lineToPrint = addParam(headerLine, lineToPrint, lineNo, "gridNumber", grids[getGridRun(mol)])

            fnAux = self._getTmpPath("tmp.mae")
            n, fnRaw = mol.poseFile.get().split('@')
            args = "-n %s %s -o %s"%(n, fnRaw, fnAux)
            self.runJob(progMaeSubset, args)

            fnOut = self._getExtraPath("pose_%d.pdb"%lineNo)
            args = "-imae %s -opdb %s"%(fnAux, fnOut)
            self.runJob(progStructConvert, args)

            headerLine, lineToPrint = addParam(headerLine, lineToPrint, lineNo, "poseFile", os.path.split(fnOut)[1])

            if lineNo==0:
                fh.write(headerLine+"\n")
            fh.write(lineToPrint+"\n")
            lineNo+=1
        fh.close()