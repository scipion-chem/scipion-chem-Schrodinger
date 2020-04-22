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
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam, FloatParam, IntParam
from pwem.protocols import EMProtocol
from schrodingerScipion import Plugin as schrodinger_plugin
from pyworkflow.utils.path import createLink, makePath

class ProtSchrodingerGlideDocking(EMProtocol):
    """Calls glide to perform a docking of a set of compounds in a pocket defined by a grid.
       It is assumed that the input library of ligands is already prepared. """
    _label = 'docking (glide)'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputGrid', PointerParam, pointerClass="SchrodingerGrid",
                       label='Grid to analyze:', allowsNull=False)
        form.addParam('inputLibrary', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Library of compounds:', allowsNull=False)
        form.addParam('dockingMethod', EnumParam, default=0, choices=['Flexible dock', 'Rigid dock',
                                                                      'Refine (do not dock)',
                                                                      'Score in place (do not dock)'],
                       label='Docking method')
        form.addParam('dockingPrecision', EnumParam, default=0, choices=['Low (HTVS)','Medium (SP)','High (XP)'],
                       label='Docking precision',
                       help='You may use a low to high strategy. HTVS takes about 2 s/ligand, SP about 10s, and XP about 10 min.')
        form.addParam('maxkeep', IntParam, default=5000, expertLevel=LEVEL_ADVANCED,
                       label='No. Poses to keep per ligand:',
                       help='Number of poses per ligand to keep in initial phase of docking.')
        form.addParam('maxkeep', IntParam, default=5000, expertLevel=LEVEL_ADVANCED,
                       label='No. Poses to keep per ligand:',
                       help='Number of poses per ligand to keep in initial phase of docking.')
        form.addParam('scoreCutoff', FloatParam, default=100.0, expertLevel=LEVEL_ADVANCED,
                       label='Score cutoff:',
                       help='Scoring window for keeping initial poses.')
        form.addParam('maxref', IntParam, default=-1, expertLevel=LEVEL_ADVANCED,
                       label='No. Poses to keep per ligand:',
                       help='Number of poses to keep per ligand for energy minimization. '
                            'If set to -1, the default value is 400, except for XP precision that is 800.')
        form.addParam('posesPerLig', IntParam, default=5, expertLevel=LEVEL_ADVANCED,
                       label='No. Poses to report per ligand:')

        form.addSection(label='Ligand sampling')
        form.addParam('sampleNinversions', BooleanParam, default=True, condition='dockingMethod==0',
                       label='Sample pyramid nitrogen inversions:')
        form.addParam('sampleRings', BooleanParam, default=True, condition='dockingMethod==0',
                       label='Sample rings:')
        form.addParam('epikPenalties', BooleanParam, default=False, label='Epik penalties:',
                      help='Apply penalties for ionization or tautomeric states calculated by Epik')
        form.addParam('skipMetalEpik', BooleanParam, default=True, label='Skip Epik metal only:',
                      help='Skip Epik-generated states of ligands that are designed for binding to metals. '
                           'This option is useful if the receptor has a metal but the ligand does not bind to it. '
                           'These states are skipped by default if the receptor does not have a metal.')
        form.addParam('expandedSampling', BooleanParam, default=False, label='Expanded sampling:',
                      help='Expand the sampling by bypassing the elimination of poses in the rough scoring stage. '
                           'Useful for fragment docking.')
        form.addParam('rewardIntraHBonds', BooleanParam, default=False, label='Reward intra H bonds:',
                      help='Reward intramolecular ligand hydrogen bonds by adding a contribution for each '
                           'intramolecular hydrogen bond to the GlideScore, and a contribution to Emodel')
        form.addParam('HbondDonorAromH', BooleanParam, default=False, label='Aromatic H as H-bond donors:',
                      help='Accept aromatic hydrogens as potential H-bond donors.')
        form.addParam('HbondDonorAromHCharge', FloatParam, default=0.0, label='Aromatic H as H-bond donors Charge:',
                      condition='HbondDonorAromH',
                      help='Partial charge cutoff for accepting aromatic hydrogens as potential H-bond donors. '
                           'The cutoff is applied to the actual (signed) charge, not the absolute value.')
        form.addParam('HbondAcceptHalo', BooleanParam, default=False, label='Halogens as H-bond acceptors:',
                      help='Accept halogens (neutral or charged, F, Cl, Br, or I) as H-bond acceptors.')
        form.addParam('HbondDonorHalo', BooleanParam, default=False, label='Halogens as H-bond donors:',
                      help='Accept the halogens (Cl, Br, I, but not F) as potential H-bond '
                           '(noncovalent interaction) donors')

        form.addParallelSection(threads=1, mpi=8)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('dockingStep')

    def dockingStep(self):
        def makeLocal(fn):
            return fn.replace(self._getPath()+'/','')

        glideProg = schrodinger_plugin.getHome('glide')

        fnGrid = self._getExtraPath("grid.zip")
        if not os.path.exists(fnGrid): # Prepared to resume
            createLink(self.inputGrid.get().getFileName(),fnGrid)

        fnLigands = self._getPath('ligands')
        if not os.path.exists(fnLigands):
            makePath(fnLigands)

        fnIn = self._getPath('job.inp')
        if not os.path.exists(fnIn): # Prepared to resume
            fhIn = open(fnIn,'w')
            fhIn.write("GRIDFILE %s\n"%makeLocal(fnGrid))

            if self.dockingMethod.get()==0:
                fhIn.write("DOCKING_METHOD confgen\n")
                fhIn.write("FLEXTORS True\n")
            elif self.dockingMethod.get()==1:
                fhIn.write("DOCKING_METHOD rigid\n")
            elif self.dockingMethod.get()==2:
                fhIn.write("DOCKING_METHOD mininplace\n")
            elif self.dockingMethod.get()==3:
                fhIn.write("DOCKING_METHOD inplace\n")

            if self.dockingPrecision.get()==0:
                fhIn.write("PRECISION HTVS\n")
            elif self.dockingPrecision.get()==1:
                fhIn.write("PRECISION SP\n")
            elif self.dockingPrecision.get()==2:
                fhIn.write("PRECISION XP\n")
                fhIn.write("WRITE_XP_DESC True\n")
                fhIn.write("POSTDOCK_NPOSE 10\n")

            fhIn.write("SAMPLE_N_INVERSIONS %s\n"%self.sampleNinversions.get())
            fhIn.write("SAMPLE_RINGS %s\n"%self.sampleRings.get())
            fhIn.write("EPIK_PENALTIES %s\n"%self.sampleNinversions.get())
            fhIn.write("SKIP_EPIK_METAL_ONLY %s\n"%self.skipMetalEpik.get())
            fhIn.write("EXPANDED_SAMPLING %s\n"%self.expandedSampling.get())
            fhIn.write("REWARD_INTRA_HBONDS %s\n"%self.rewardIntraHBonds.get())
            fhIn.write("HBOND_DONOR_AROMH %s\n"%self.HbondDonorAromH.get())
            if self.HbondDonorAromH.get():
                fhIn.write("HBOND_DONOR_AROMH_CHARGE %f\n" % self.HbondDonorAromHCharge.get())
            fhIn.write("HBOND_ACCEP_HALO %s\n"%self.HbondAcceptHalo.get())
            fhIn.write("HBOND_DONOR_HALO %s\n"%self.HbondDonorHalo.get())

            fhIn.write("MAXKEEP %d\n"%self.maxkeep.get())
            fhIn.write("SCORING_CUTOFF %f\n"%self.scoreCutoff.get())
            if self.maxref.get()>0:
                fhIn.write("MAXREF %d\n" % self.maxref.get())
            else:
                if self.dockingPrecision.get()==2:
                    fhIn.write("MAXREF %d\n" % 800)
                else:
                    fhIn.write("MAXREF %d\n" % 400)
            fhIn.write("WRITE_CSV True\n")
            fhIn.write("POSES_PER_LIG %d\n"%self.posesPerLig.get())

            fhIn.write("LIGANDFILE ")
            i=0
            for small in self.inputLibrary.get():
                fnSmall = small.getFileName()
                fnSmallLocal = os.path.join(fnLigands,os.path.split(fnSmall)[1])
                createLink(fnSmall,fnSmallLocal)
                if i!=0:
                    fhIn.write(",")
                fhIn.write(makeLocal(fnSmallLocal))
                i+=1
            fhIn.write("\n")

            fhIn.close()

        args = "-WAIT -NJOBS %d -RESTART -LOCAL job.inp"%(self.numberOfMpi.get())

        self.runJob(glideProg, args, cwd=self._getPath())

    def _validate(self):
        errors = []
        if self.HbondAcceptHalo.get() and self.HbondDonorHalo.get():
            errors.append('Halogens cannot be simultaneously acceptors and donors of H-bonds')
        return errors