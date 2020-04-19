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

from pyworkflow.protocol.params import PointerParam, EnumParam
from pwem.objects.data import AtomStruct
from pwem.protocols import EMProtocol
from schrodinger import Plugin
from schrodinger.objects import SchrodingerAtomStruct
from bioinformatics.objects import SetOfSmallMolecules, SmallMolecule

class ProtSchrodingerConvert(EMProtocol):
    """Convert a set of input ligands or a receptor structure to a specific file format"""
    _label = 'convert'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputType', EnumParam, default=0,
                       choices=["Small molecules",'Target structure'],
                       label='Input type')
        form.addParam('inputSmallMols', PointerParam, pointerClass="SetOfSmallMolecules", condition='inputType==0',
                      label='Input set:', allowsNull=False)
        form.addParam('outputFormatSmall', EnumParam, default=2, condition='inputType==0',
                       choices=["Maestro",'PDB',"Sybyl Mol2","Smiles",'Cif','V2000 SD','CSV (Smiles with properties)'],
                       label='Output format')
        form.addParam('inputStructure', PointerParam, pointerClass="SchrodingerAtomStruct, AtomStruct",
                      condition='inputType==1', label='Input structure:', allowsNull=False)
        form.addParam('outputFormatTarget', EnumParam, default=2, condition='inputType==1',
                       choices=["Maestro",'PDB',"Sybyl Mol2"],
                       label='Output format')

        # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertStep')

    def convertStep(self):
        def inputArg(fn):
            if fn.endswith('.mae') or fn.endswith('.maegz'):
                args="-imae "
            elif fn.endswith('.cif'):
                args="-icif"
            elif fn.endswith('.sd'):
                args="-isd"
            elif fnSmall.endswith('.pdb'):
                args="-ipdb"
            elif fnSmall.endswith('.mol2'):
                args="-imol2"
            elif fn.endswith('.smi'):
                args="-ismi"
            elif fnSmall.endswith('.csv'):
                args="-icsv"
            return args

        def outputArg(fnRoot, format):
            if format==0:
                fnOut=self._getExtraPath(fnRoot+".maegz")
                args=" -omae %s"%fnOut
            elif format == 1:
                fnOut = self._getExtraPath(fnRoot + ".pdb")
                args = " -opdb %s" % fnOut
            elif format == 2:
                fnOut = self._getExtraPath(fnRoot + ".mol2")
                args = " -omol2 %s" % fnOut
            elif format==3:
                fnOut=self._getExtraPath(fnRoot + ".smi")
                args= " -osmi %s"%fnOut
            elif format==4:
                fnOut=self._getExtraPath(fnRoot + ".cif")
                args= " -ocif %s"%fnOut
            elif format==5:
                fnOut=self._getExtraPath(fnRoot + ".sd")
                args= " -osd %s"%fnOut
            elif format==6:
                fnOut=self._getExtraPath(fnRoot + ".csv")
                args= " -ocsv %s"%fnOut
            return fnOut,args

        progStructConvert=Plugin.getHome('utilities/structconvert')

        if self.inputType==0:
            outputSmallMolecules = SetOfSmallMolecules().create(path=self._getPath(),suffix='SmallMols')

            for mol in self.inputSmallMols.get():
                fnSmall = mol.smallMoleculeFile.get()
                fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]

                args=inputArg(fnSmall)+" %s"%fnSmall
                fnOut, argout = outputArg(fnRoot, self.outputFormatSmall.get())
                args+=argout

                self.runJob(progStructConvert, args)
                smallMolecule = SmallMolecule(smallMolFilename=fnOut)
                outputSmallMolecules.append(smallMolecule)

            if len(outputSmallMolecules)>0:
                self._defineOutputs(outputSmallMols=outputSmallMolecules)
                self._defineSourceRelation(self.inputSmallMols, outputSmallMolecules)
        else:
            fnStructure = self.inputStructure.get().getFileName()
            args = inputArg(fnStructure) + " %s" % fnStructure
            fnRoot = os.path.splitext(os.path.split(fnStructure)[1])[0]
            fnOut, argout = outputArg(fnRoot, self.outputFormatTarget.get())
            args += argout
            self.runJob(progStructConvert, args)

            if fnOut.endswith('.maegz'):
                target = SchrodingerAtomStruct(filename=fnOut)
            else:
                target = AtomStruct(filename=fnOut)
            self._defineOutputs(outputStructure=target)
            self._defineSourceRelation(self.inputStructure, target)

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