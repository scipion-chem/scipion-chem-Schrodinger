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
import os, time

from pyworkflow.protocol.params import PointerParam, EnumParam, STEPS_PARALLEL
from pwem.objects.data import AtomStruct
from pwem.protocols import EMProtocol
from schrodingerScipion import Plugin
from schrodingerScipion.objects import SchrodingerAtomStruct
from schrodingerScipion.utils.utils import putMol2Title
from pwchem.objects import SetOfSmallMolecules, SmallMolecule

SMALLMOL, TARGET = 0, 1

class ProtSchrodingerConvert(EMProtocol):
    """Convert a set of input ligands or a receptor structure to a specific file format"""
    _label = 'convert'
    _program = ""
    _outFormats = ['maegz', 'pdb', 'mol2', 'smi', 'cif', 'sd', 'csv']

    saving = False

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputType', EnumParam, default=0,
                       choices=["Small molecules", 'Target structure'],
                       label='Input type')
        form.addParam('inputSmallMols', PointerParam, pointerClass="SetOfSmallMolecules",
                      condition='inputType=={}'.format(SMALLMOL),
                      label='Input set:', allowsNull=False)
        form.addParam('outputFormatSmall', EnumParam, default=2, condition='inputType=={}'.format(SMALLMOL),
                       choices=["Maestro", 'PDB', "Sybyl Mol2", "Smiles", 'Cif',
                                'V2000 SD', 'CSV (Smiles with properties)'],
                       label='Output format')
        form.addParam('inputStructure', PointerParam, pointerClass="SchrodingerAtomStruct, AtomStruct",
                      condition='inputType=={}'.format(TARGET), label='Input structure:', allowsNull=False)
        form.addParam('outputFormatTarget', EnumParam, default=2, condition='inputType=={}'.format(TARGET),
                       choices=["Maestro", 'PDB', "Sybyl Mol2"],
                       label='Output format')

        form.addParallelSection(threads=4, mpi=1)

        # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        convSteps = []
        if self.inputType.get() == SMALLMOL:
            self.outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')
            for mol in self.inputSmallMols.get():
                cStep = self._insertFunctionStep('convertMolStep', mol.clone(), prerequisites=[])
                convSteps.append(cStep)

        elif self.inputType.get() == TARGET:
            cStep = self._insertFunctionStep('convertTargetStep', prerequisites=[])
            convSteps.append(cStep)

        self._insertFunctionStep('createOutputStep', prerequisites=convSteps)

    def convertMolStep(self, mol):
        progStructConvert=Plugin.getHome('utilities/structconvert')

        fnSmall = mol.smallMoleculeFile.get()
        fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]
        fnOut = self._getExtraPath("{}.{}".format(fnRoot, self._outFormats[self.outputFormatSmall.get()]))

        args = '{} {}'.format(fnSmall, fnOut)

        self.runJob(progStructConvert, args)
        if self.outputFormatSmall.get()==2:
            putMol2Title(fnOut)

        self.saveMolecule(fnOut, self.outputSmallMolecules)

    def convertTargetStep(self):
        progStructConvert = Plugin.getHome('utilities/structconvert')
        fnStructure = self.inputStructure.get().getFileName()
        args = fnStructure
        fnRoot = os.path.splitext(os.path.split(fnStructure)[1])[0]
        fnOut = self._getExtraPath("{}.{}".format(fnRoot, self._outFormats[self.outputFormatTarget.get()]))
        args += fnOut
        self.runJob(progStructConvert, args)

        if fnOut.endswith('.maegz'):
            self.target = SchrodingerAtomStruct(filename=fnOut)
        else:
            self.target = AtomStruct(filename=fnOut)

    def createOutputStep(self):
        if self.inputType.get() == 0:
            self._defineOutputs(outputSmallMols=self.outputSmallMolecules)
        else:
            self._defineOutputs(outputStructure=self.target)

    def saveMolecule(self, molFn, molSet):
        while self.saving:
            time.sleep(0.2)
        self.saving = True
        smallMolecule = SmallMolecule(smallMolFilename=molFn)
        molSet.append(smallMolecule.clone())
        self.saving = False


    def _summary(self):
        summary=[]
        if self.inputType.get()==0:
            if self.outputFormatSmall.get() == 0:
                summary.append('Converted to Maestro')
            elif self.outputFormatTarget.get() == 1:
                summary.append('Converted to PDB')
            elif self.outputFormatSmall.get() == 2:
                summary.append('Converted to Sybyl Mol2')
            elif self.outputFormatSmall.get() == 3:
                summary.append('Converted to Smiles')
            elif self.outputFormatSmall.get() == 4:
                summary.append('Converted to Cif')
            elif self.outputFormatSmall.get() == 5:
                summary.append('Converted to V2000 SD')
            elif self.outputFormatSmall.get() == 6:
                summary.append('Converted to CSV (Smiles with properties)')
        else:
            if self.outputFormatTarget.get() == 0:
                summary.append('Converted to Maestro')
            elif self.outputFormatTarget.get() == 1:
                summary.append('Converted to PDB')
            elif self.outputFormatTarget.get() == 2:
                summary.append('Converted to Sybyl Mol2')
        return summary