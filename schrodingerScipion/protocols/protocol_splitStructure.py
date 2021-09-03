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

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam, StringParam
from pwem.protocols import EMProtocol
from schrodingerScipion import Plugin
from schrodingerScipion.objects import SchrodingerAtomStruct
from pwchem.objects import SmallMolecule

class ProtSchrodingerSplitStructure(EMProtocol):
    """Split a structure into different pieces"""
    _label = 'split structure'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="SchrodingerAtomStruct",
                      label='Input structure:', allowsNull=False)
        form.addParam('splitMode', EnumParam, default=0,
                       choices=["Chain","Molecule","Residue","Ligand","PDB"],
                       label='Split mode',
                       help='Mode=pdb will split the structure into 1) receptor, 2) each individual ligand, '
                            '3) all non-metal ions and cofactors, 4) all waters. By default, mode=chain. '
                            'Use mode=ligand to remove ligands')
        form.addParam('mergeLigands', BooleanParam, default=False, condition="splitMode==0",
                      label='Merge ligands with the closest chain')
        form.addParam('mergeWaters', BooleanParam, default=False, condition="splitMode==0",
                      label='Merge waters with the closest chain')
        form.addParam('keepProperties', BooleanParam, default=False,
                      label='Keep properties')
        form.addParam('groupWaters', BooleanParam, default=False, condition="splitMode!=4",
                      label='Merge waters with the closest chain')
        form.addParam('splitAll', BooleanParam, default=False,
                      label='Split cofactors and metals into different structures')
        form.addParam('ligandASL', StringParam, default="", expertLevel=LEVEL_ADVANCED,
                      label='ASL used to define ligand structures',
                      help='For help on ASL (Atom Specification Language), see Chap. 3 of '
                           'http://shaker.umh.es/computing/Schrodinger_suites/maestro_command_reference.pdf')
        form.addParam('cofactorASL', StringParam, default="", expertLevel=LEVEL_ADVANCED,
                      label='ASL used to define cofactor structures',
                      help='For help on ASL (Atom Specification Language), see Chap. 3 of '
                           'http://shaker.umh.es/computing/Schrodinger_suites/maestro_command_reference.pdf')
        form.addParam('positiveASL', StringParam, default="", expertLevel=LEVEL_ADVANCED, condition="splitMode==4",
                      label='ASL used to define non-metal positive ions',
                      help='For help on ASL (Atom Specification Language), see Chap. 3 of '
                           'http://shaker.umh.es/computing/Schrodinger_suites/maestro_command_reference.pdf')
        form.addParam('negativeASL', StringParam, default="", expertLevel=LEVEL_ADVANCED, condition="splitMode==4",
                      label='ASL used to define negative ions',
                      help='For help on ASL (Atom Specification Language), see Chap. 3 of '
                           'http://shaker.umh.es/computing/Schrodinger_suites/maestro_command_reference.pdf')

        # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('splitStep')

    def splitStep(self):
        args=Plugin.getMMshareDir('python/common/split_structure.py')
        if self.splitMode.get()==0:
            args+=" -m chain"
            if self.mergeLigands.get():
                args+=" -merge_ligands_with_chain"
            if self.mergeWaters.get():
                args+=" -merge_waters_with_chain"
        elif self.splitMode.get()==1:
            args+=" -m molecule"
        elif self.splitMode.get()==2:
            args+=" -m residue"
        elif self.splitMode.get()==3:
            args+=" -m ligand"
        elif self.splitMode.get()==4:
            args+=" -m pdb"
            if self.positiveASL.get()!="":
                args+=' -positive_ion_asl "%s"'%self.positiveASL.get()
            if self.negativeASL.get()!="":
                args+=' -negative_ion_asl "%s"'%self.negativeASL.get()
        if self.keepProperties.get():
            args+=" --keep_properties"
        if self.splitMode.get()!=4 and self.groupWaters.get():
            args+=" -groupwaters"
        if self.splitAll.get():
            args+=" -splitall"
        if self.ligandASL.get()!="":
            args+=' -ligand_asl "%s"'%self.ligandASL.get()
        if self.cofactorASL.get()!="":
            args+=' -cofactor_asl "%s"'%self.cofactorASL.get()
        args+=" -many_files %s %s"%(self.inputStructure.get().getFileName(),
                                    self._getExtraPath('output.maegz'))

        self.runJob(Plugin.getHome('run'),args)
        def getNumber(fn,suffix):
            fnBase = os.path.splitext(os.path.split(fn)[1])[0]
            tokens = fnBase.split(suffix)
            return tokens[1]

        for fn in glob.glob(self._getExtraPath("output*")):
            if "_receptor" in fn:
                target = SchrodingerAtomStruct(filename=fn)
                number = getNumber(fn,"_receptor")
                outputDict = {'outputStructure%s' % number: target}
                self._defineOutputs(**outputDict)
                self._defineSourceRelation(self.inputStructure, target)
            elif "_ligand" in fn:
                ligand = SmallMolecule(smallMolFilename=fn)
                number = getNumber(fn, "_ligand")
                outputDict = {'outputLigand%s' % number: ligand}
                self._defineOutputs(**outputDict)
                self._defineSourceRelation(self.inputStructure, ligand)
            elif "_cof_ion" in fn:
                cof = SchrodingerAtomStruct(filename=fn)
                number = getNumber(fn, "_cof_ion")
                outputDict = {'outputCofactors%s' % number: cof}
                self._defineOutputs(**outputDict)
                self._defineSourceRelation(self.inputStructure, cof)
            else:
                print("Scipion: I don't know how to handle %s"%fn)

    def _summary(self):
        summary=[]
        return summary