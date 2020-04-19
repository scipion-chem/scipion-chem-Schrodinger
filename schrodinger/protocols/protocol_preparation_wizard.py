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
from pyworkflow.protocol.params import PointerParam, BooleanParam, FloatParam, IntParam, EnumParam
from pyworkflow.utils.path import createLink
from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct
from pwem.convert.atom_struct import AtomicStructHandler
from schrodinger import Plugin
from schrodinger.objects import SchrodingerAtomStruct

class ProtSchrodingerPrepWizard(EMProtocol):
    """Calls the preparation wizard"""
    _label = 'target preparation (prepwizard)'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct, SchrodingerAtomStruct",
                       label='Atomic Structure:', allowsNull=False)
        form.addSection(label='Stage 1')
        form.addParam('stage1', BooleanParam, default=False, label='Stage 1 preparation:')
        form.addParam('fillSideChains', BooleanParam, default=True, condition='stage1',
                      label='Fill side chains:',
                      help='Fill missing side chains Prime')
        form.addParam('fillLoops', BooleanParam, default=True, condition='stage1',
                      label='Fill loops:',
                      help='Fill missing loops with Prime')
        form.addParam('disulfides', BooleanParam, default=True, condition='stage1',
                      label='Create bonds to proximal Sulfurs:',
                      help='Delete hydrogens as needed')
        form.addParam('mse', BooleanParam, default=True, condition='stage1',
                      label='Convert Selenomethionine residues to Methionines:')
        form.addParam('hydrogens', EnumParam, default=0, choices=["Don't add hydrogens", "Delete and re-add hydrogens"],
                      condition='stage1',
                      label='Hydrogen treatment:')
        form.addParam('glycosylation', BooleanParam, default=False, condition='stage1',
                      label='Glycosylation:',
                      help='Create bonds to N-linked and O-linked sugars (delete hydrogens as needed)')
        form.addParam('palmitoylation', BooleanParam, default=False, condition='stage1',
                      label='Palmitoylation:',
                      help='Create palmitoylation bonds even if not included in the CONNECT records ( delete hydrogens as needed)')
        form.addParam('captermini', BooleanParam, default=True, condition='stage1',
                      label='Add cap termini:',
                      help='Add ACE and NME termini')
        form.addParam('keepFarWat', BooleanParam, default=False, condition='stage1',
                      label='Keep far waters:',
                      help="Don't delete waters far from het groups")
        form.addParam('watdist', FloatParam, default=5, condition='stage1 and keepFarWat',
                      label='Water distance cutoff (A)',
                      help='Distance threshold for far waters')
        form.addParam('treatMetals', BooleanParam, default=True, condition='stage1',
                      label='Treat metals:',
                      help="Adding zero-order bonds, etc")

        form.addSection(label='Stage 2. Protonation (Protassign)')
        form.addParam('stage2', BooleanParam, default=False, label='Stage 2 preparation:')
        form.addParam('sampleWaters', BooleanParam, default=True, condition='stage2',
                      label='Sample waters:')
        form.addParam('xtal', BooleanParam, default=False, condition='stage2',
                      label='Use crystal symmetry:')
        form.addParam('propKa', BooleanParam, default=False, condition='stage2',
                      label='Adjust protonation to a specific pH:')
        form.addParam('propKapH', FloatParam, default=7, condition='stage2 and propKa',
                      label='pH:')
        form.addParam('minadjh', BooleanParam, default=True, condition='stage2',
                      label='Minimize all adjustable hydrogens:',
                      help='Titratable hydrogens, water hydrogens, hydroxyls, thiols, and ASN/GLN carboxamide hydrogens')

        form.addSection(label='Stage 3. Restrained minimization (Impref)')
        form.addParam('stage3', BooleanParam, default=False, label='Stage 3 preparation:')
        form.addParam('rmsdD', FloatParam, default=0.3, condition='stage3',
                      label='RMSD cutoff:')
        form.addParam('fix', BooleanParam, default=False, condition='stage3',
                      label='Fix heavy atoms:', help='Minimize hydrogens only')
        form.addParam('force', EnumParam, default=0, choices=["2005","3"],
                      condition='stage3', label='OPLS force field:')

        form.addSection(label='Stage 4. Ionization and tautomeric states of HET groups (Epik)')
        form.addParam('stage4', BooleanParam, default=False, label='Stage 4 preparation:')
        form.addParam('ms', BooleanParam, default=False, condition='stage4',
                      label='Apply Max. States constrain:')
        form.addParam('msN', IntParam, default=1, condition='stage4 and ms',
                      label='Max. Number of states:')
        form.addParam('epikPh', FloatParam, default=7, condition='stage4',
                      label='pH:')
        form.addParam('epikPht', FloatParam, default=2, condition='stage4',
                      label='pH Range:')


    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('preparationStep')
        self._insertFunctionStep('createOutput')

    def preparationStep(self):
        prog=Plugin.getHome('utilities/prepwizard')

        if isinstance(self.inputStructure.get(),AtomStruct):
            fnIn = self._getExtraPath("atomStructIn.pdb")
            aStruct1 = AtomicStructHandler(self.inputStructure.get().getFileName())
            aStruct1.write(fnIn)
            fnIn='extra/atomStructIn.pdb'
        else:
            fnIn = self._getExtraPath("atomStructIn.mae")
            createLink(self.inputStructure.get().getFileName(),fnIn)
            fnIn='extra/atomStructIn.mae'

        args='-WAIT'
        if self.stage1.get():
            if self.fillSideChains.get():
                args+=' -fillsidechains'
            if self.fillLoops.get():
                args+=' -fillloops'
            if self.disulfides.get():
                args+=' -disulfides'
            if self.mse.get():
                args+=' -mse'
            if self.hydrogens.get()==0:
                args+=" -nohtreat"
            elif self.hydrogens.get()==1:
                args+=" -rehtreat"
            if self.glycosylation.get():
                args+=" -glycosylation"
            if self.palmitoylation.get():
                args+=" -palmitoylation"
            if self.captermini.get():
                args+=" -captermini"
            if self.keepFarWat.get():
                args+=" -keepfarwat -watdist %f"%self.watdist.get()
            if not self.treatMetals.get():
                args+=" -nometaltreat"

        if self.stage2.get():
            if self.sampleWaters.get():
                args+=" -samplewater"
            if self.xtal.get():
                args+=" -xtal"
            if self.propKa.get():
                args+=" -propka_pH %f"%self.propKapH.get()
            else:
                args+=" -nopropka"
            if self.minadjh.get():
                args+=" -minimize_adj_h"
        else:
            args+=" -noprotassign"

        if self.stage3.get():
            args+=" -rmsd %f"%self.rmsdD.get()
            if self.fix.get():
                args+=" -fix"
            if self.force.get()==0:
                args+=" -f 2005"
            else:
                args+=" -f 3"
        else:
            args+=" -noimpref"

        if self.stage4.get():
            if self.ms.get():
                args+=" -ms %d"%self.msN.get()
            args+=" -epik_pH %f"%self.epikPh.get()
            args+=" -epik_pHt %f"%self.epikPht.get()
        else:
            args+=" -noepik"

        args+=' %s atomStructOut.maegz'%fnIn
        self.runJob(prog,args,cwd=self._getPath())

    def createOutput(self):
        fnMae = self._getPath('atomStructOut.maegz')
        if os.path.exists(fnMae):
            maeFile=SchrodingerAtomStruct()
            maeFile.setFileName(fnMae)

            self._defineOutputs(outputStructure=maeFile)
            self._defineSourceRelation(self.inputStructure, maeFile)

    def _citations(self):
        return ['Sastry2013']
