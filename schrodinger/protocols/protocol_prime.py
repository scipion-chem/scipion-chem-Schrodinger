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

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, EnumParam, FloatParam
from pyworkflow.utils.path import createLink
from pwem.protocols import EMProtocol
from schrodinger import Plugin
from schrodinger.objects import SchrodingerAtomStruct

class ProtSchrodingerPrime(EMProtocol):
    """Schrodinger's prime is a structure prediction program """
    _label = 'prime'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="SchrodingerAtomStruct",
                       label='Atomic Structure:', allowsNull=False)
        form.addParam('operation', EnumParam, choices=['Side chain prediction','Minimization of all hydrogens',
                                                       'Loop prediction'],
                      label='Operation', default=0)
        form.addParam('residueList', StringParam, condition='operation==0',
                      label='List of residues:',
                      help='Example: Residues 110 of chain A and 45 of chain B -> A:110, B:45')
        form.addParam('residueFirst', StringParam, condition='operation==2',
                      label='First residue of the loop:',
                      help='Example: Residue 110 of chain A -> A:110')
        form.addParam('residueLast', StringParam, condition='operation==2',
                      label='Last residue of the loop:',
                      help='Example: Residue 115 of chain A -> A:115')
        form.addParam('resSphere', FloatParam, condition='operation==2', expertLevel=LEVEL_ADVANCED,
                      label='Residue sphere (A):', default=7.5)
        form.addParam('minOverlap', FloatParam, condition='operation==2', expertLevel=LEVEL_ADVANCED,
                      label='Minimum overlap:', default=0.7)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('primeStep')
        self._insertFunctionStep('createOutput')

    def primeStep(self):
        prog=Plugin.getHome('prime')

        fnInputStructure = self.inputStructure.get().getFileName()
        if fnInputStructure.endswith('.mae'):
            fnIn = self._getPath("atomStructIn.mae")
        elif fnInputStructure.endswith('.maegz'):
            fnIn = self._getPath("atomStructIn.maegz")

        createLink(fnInputStructure,fnIn)
        fnIn=os.path.split(fnIn)[1]

        fhJob = open(self._getPath('job.inp'),'w')
        fhJob.write("STRUCT_FILE %s\n" % fnIn)
        if self.operation.get()==0:
            # Side chain prediction
            fhJob.write("PRIME_TYPE  SIDE_PRED\n")
            fhJob.write("SELECT  pick\n")
            i=0
            for residue in self.residueList.get().split(','):
                fhJob.write("RESIDUE_%d %s\n"%(i,residue.strip()))
                i+=1
        elif self.operation.get()==1:
            # Minimization of all hydrogens
            fhJob.write("PRIME_TYPE  REAL_MIN\n")
            fhJob.write("SELECT  asl = (atom.ele H)\n")
        elif self.operation.get()==2:
            # Loop prediction
            fhJob.write("PRIME_TYPE  LOOP_BLD\n")
            fhJob.write("LOOP_0_RES_0 %s\n"%self.residueFirst.get())
            fhJob.write("LOOP_0_RES_1 %s\n"%self.residueLast.get())
            fhJob.write("RES_SPHERE %f\n"%self.resSphere.get())
            fhJob.write("MIN_OVERLAP %f\n"%self.minOverlap.get())

        fhJob.write("USE_CRYSTAL_SYMMETRY no\n")
        fhJob.write("USE_RANDOM_SEED no\n")
        fhJob.write("SEED 0\n")
        fhJob.write("EXT_DIEL 80.00\n")
        fhJob.write("USE_MEMBRANE no\n")
        fhJob.close()
        args='-WAIT -LOCAL job.inp'
        self.runJob(prog,args,cwd=self._getPath())

    def createOutput(self):
        fnMae = self._getPath('job-out.maegz')
        if os.path.exists(fnMae):
            maeFile=SchrodingerAtomStruct()
            maeFile.setFileName(fnMae)

            self._defineOutputs(outputStructure=maeFile)
            self._defineSourceRelation(self.inputStructure, maeFile)

    def _citations(self):
        return ['Jacobson2002', 'Jacobson2004']
