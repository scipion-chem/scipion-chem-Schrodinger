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

from pyworkflow.protocol.params import PointerParam, STEPS_PARALLEL, FloatParam, BooleanParam, IntParam
from pwem.protocols import EMProtocol
from schrodingerScipion import Plugin
from pwchem.utils import mergePDBs

progStructConvert=Plugin.getHome('utilities/structconvert')
progMaeSubset=Plugin.getHome('utilities/maesubset')

class ProtSchrodingerPrepareForWeb(EMProtocol):
    """Prepare a set of small molecules for the web"""
    _label = 'prepare for web'
    _program = ""

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSmallMols', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Set of small molecules: ', allowsNull=False)
        form.addParam('addTarget', BooleanParam, default=False,
                      label='Add protein structure: ',
                      help='Add the target structure to the output pdbs')
        form.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                      label='Optional atomic structure: ', condition='addTarget', allowsNull=True,
                      help='Target structure to include in the output pdbs')
        form.addParallelSection(threads=4, mpi=1)

        # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self.csvStr, prepSteps = '', []
        for i, mol in enumerate(self.inputSmallMols.get()):
            pStep = self._insertFunctionStep('prepareStep', mol.clone(), i, prerequisites=[])
            prepSteps.append(pStep)

        self._insertFunctionStep('createCSVStep', prerequisites=prepSteps)

    def createCSVStep(self):
        with open(self._getExtraPath("smallMolecules.csv"), 'w') as fh:
            fh.write(self.csvStr)

    def prepareStep(self, mol, lineNo):
        if self.addTarget:
            os.mkdir(self._getExtraPath('withTarget'))

        headerLine = ""
        fnSmall = mol.smallMoleculeFile.get()
        fnMol = os.path.split(fnSmall)[1]
        fnRoot = os.path.splitext(fnMol)[0]
        if "-" in fnRoot:
            fnRoot = fnRoot.split('-')[0]

        lineToPrint = ""
        headerLine, lineToPrint = self.addMolParam(mol, headerLine, lineToPrint, lineNo, "smallMolecule", fnRoot)
        if fnRoot.startswith("ZINC"):
            headerLine, lineToPrint = self.addMolParam(mol, headerLine, lineToPrint, lineNo, "url",
                                                    "https://zinc15.docking.org/substances/%s"%fnRoot)
        headerLine, lineToPrint = self.addMolParam(mol, headerLine, lineToPrint, lineNo, "bindingSiteScore")
        headerLine, lineToPrint = self.addMolParam(mol, headerLine, lineToPrint, lineNo, "dockingScore")
        headerLine, lineToPrint = self.addMolParam(mol, headerLine, lineToPrint, lineNo, "ligandEfficiency")

        headerLine, lineToPrint = self.addMolParam(mol, headerLine, lineToPrint, lineNo, "gridNumber",
                                                self.getGridId(mol))

        fnOut = self._getExtraPath("pose_%d.pdb"%lineNo)
        if not os.path.exists(fnOut):
            fnAux = self._getTmpPath("tmp_%d.mae"%lineNo)
            n, fnRaw = mol.poseFile.get().split('@')
            args = "-n %s %s -o %s"%(n, fnRaw, fnAux)
            self.runJob(progMaeSubset, args)

            args = "%s %s"%(fnAux, fnOut)
            self.runJob(progStructConvert, args)

        if self.addTarget:
            mergePDBs(self.inputAtomStruct.get().getFileName(), fnOut,
                      self._getExtraPath('withTarget/targetPose_%d.pdb' % lineNo))

        headerLine, lineToPrint = self.addMolParam(mol, headerLine, lineToPrint, lineNo, "poseFile",
                                                os.path.split(fnOut)[1])
        if lineNo==0:
            self.csvStr += headerLine+"\n"
        self.csvStr += lineToPrint+"\n"

    def addMolParam(self, mol, headerLine, lineToPrint, lineNo, param, value=None):
        if value is None:
            if hasattr(mol, param):
                value = getattr(mol, param).get()
        if value is not None:
            if lineNo == 0:
                if headerLine != "":
                    headerLine += "; "
                headerLine += param
            if lineToPrint != "":
                lineToPrint += "; "
            lineToPrint += str(value)
        return headerLine, lineToPrint

    def getGridId(self, mol):
        fn = mol.poseFile.get()
        return fn.split('/')[-2].split('grid_')[1]