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
import os, shutil, threading, subprocess

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam, FloatParam, IntParam, STEPS_PARALLEL
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.utils.path import makePath
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from schrodingerScipion import Plugin as schrodinger_plugin
from schrodingerScipion.utils.utils import putMol2Title, sortDockingResults
from schrodingerScipion.objects import SchrodingerPoses

glideProg = schrodinger_plugin.getHome('glide')
structConvertProg = schrodinger_plugin.getHome('utilities/structconvert')
structCatProg = schrodinger_plugin.getHome('utilities/structcat')
propListerProg = schrodinger_plugin.getHome('utilities/proplister')
maeSubsetProg = schrodinger_plugin.getHome('utilities/maesubset')

class ProtSchrodingerGlideDocking(EMProtocol):
    """Calls glide to perform a docking of a set of compounds in a pocket defined by a grid.
       It is assumed that the input library of ligands is already prepared.

       The dockinsScore is the Glide Docking Score and it is measured in kcal/mol"""
    _label = 'docking (glide)'
    _program = ""

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputGridSet', PointerParam, pointerClass="SetOfSchrodingerGrids",
                       label='Grid to analyze:', allowsNull=False)
        form.addParam('inputLibrary', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Library of compounds:', allowsNull=False)
        form.addParam('mergeOutput', BooleanParam, default=False,
                      label='Merge output from different pockets: ')
        form.addParam('doConvertOutput', BooleanParam, default=False,
                      label='Convert output from maestro format: ')
        form.addParam('convertType', EnumParam, condition='doConvertOutput',
                      choices=['pdb', 'mol2', 'sdf'], default=0, display=EnumParam.DISPLAY_HLIST,
                      label='Convert to: ',
                      help='Convert output molecules to this format')
        group = form.addGroup('Docking')
        group.addParam('dockingMethod', EnumParam, default=0, choices=['Flexible dock (confgen)', 'Rigid dock (rigid)',
                                                                      'Refine (do not dock, mininplace)',
                                                                      'Score in place (do not dock, inplace)'],
                       label='Docking method')
        group.addParam('dockingPrecision', EnumParam, default=0, choices=['Low (HTVS)','Medium (SP)','High (XP)'],
                       label='Docking precision',
                       help='You may use a low to high strategy. HTVS takes about 2 s/ligand, SP about 10s, and XP about 10 min.')
        group.addParam('maxkeep', IntParam, default=5000, expertLevel=LEVEL_ADVANCED,
                       label='No. Poses to keep per ligand:',
                       help='Number of poses per ligand to keep in initial phase of docking.')
        group.addParam('scoreCutoff', FloatParam, default=100.0, expertLevel=LEVEL_ADVANCED,
                       label='Score cutoff:',
                       help='Scoring window for keeping initial poses.')
        group.addParam('maxref', IntParam, default=-1, expertLevel=LEVEL_ADVANCED,
                       label='No. Poses to keep per ligand:',
                       help='Number of poses to keep per ligand for energy minimization. '
                            'If set to -1, the default value is 400, except for XP precision that is 800.')
        group.addParam('posesPerLig', IntParam, default=5, expertLevel=LEVEL_ADVANCED,
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

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        cStep = self._insertFunctionStep('createLigandsFileStep', prerequisites=[])
        c2Step = self._insertFunctionStep('combineLigandsFilesStep', prerequisites=[cStep])
        dockSteps = []
        for grid in self.inputGridSet.get():
            dStep = self._insertFunctionStep('dockingStep', grid.clone(), prerequisites=[c2Step])
            dockSteps.append(dStep)
        self._insertFunctionStep('createOutput', prerequisites=dockSteps)

    def createLigandsFileStep(self):
        nt = self.numberOfThreads.get()
        ligSet = self.inputLibrary.get()
        nLigs = len(ligSet) // nt
        it, curLigSet = 0, []
        threads = []
        for lig in ligSet:
            curLigSet.append(lig.clone())
            if len(curLigSet) == nLigs and it < nt -1:
                t = threading.Thread(target=self.createLigandsFile, args=(curLigSet.copy(), it,),
                                     daemon=False)
                threads.append(t)
                t.start()
                curLigSet, it = [], it + 1

        if len(curLigSet) > 0:
            t = threading.Thread(target=self.createLigandsFile, args=(curLigSet.copy(), it,),
                                 daemon=False)
            threads.append(t)
            t.start()

        for t in threads:
            t.join()

    def combineLigandsFilesStep(self):
        nt = self.numberOfThreads.get()
        allLigandsFile = self._getTmpPath(self.getAllLigandsFile())
        with open(allLigandsFile, 'w') as f:
            for it in range(nt):
                with open(self._getTmpPath(self.getAllLigandsFile(suffix=it))) as fLig:
                    f.write(fLig.read())


    def createLigandsFile(self, ligSet, it):
        curAllLigandsFile = self._getTmpPath(self.getAllLigandsFile(suffix=it))
        with open(curAllLigandsFile, 'w') as fh:
            for small in ligSet:
                fnSmall = small.getFileName()
                if not fnSmall.endswith('.mol2'):
                    fnSmall = self.convert2Mol2(fnSmall, it)

                with open(fnSmall) as fhLigand:
                    fh.write(fhLigand.read())

    def dockingStep(self, grid):
        gridId = grid.getObjId()
        gridDir = 'grid_{}/'.format(gridId)

        makePath(self._getPath(gridDir))
        fnGrid = self._getPath(gridDir + "grid.zip")
        if not os.path.exists(fnGrid): # Prepared to resume
            shutil.copy(grid.getFileName(), fnGrid)

        fnIn = self._getPath(gridDir + 'job_{}.inp'.format(gridId))
        if not os.path.exists(fnIn): # Prepared to resume
            with open(fnIn, 'w') as fhIn:
                fhIn.write("GRIDFILE %s\n" % ("grid.zip"))

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
                fhIn.write("POSES_PER_LIG %d\n"%self.posesPerLig.get())

                fhIn.write("LIGANDFILE ../tmp/{}\n".format(self.getAllLigandsFile()))

        args = "-WAIT -RESTART -LOCAL job_{}.inp".format(gridId)
        self.runJob(glideProg, args, cwd=self._getPath(gridDir))

        self.runJob(propListerProg,
                    '-p "title" -p "docking score" -p "glide ligand efficiency" -p "glide ligand efficiency sa" -p "glide ligand efficiency ln" -c -o %s %s'%\
                    ("job_{}_pv.csv".format(gridId), "job_{}_pv.maegz".format(gridId)),
                    cwd=self._getPath(gridDir))

    def createOutput(self):
        smallDict = {}
        for small in self.inputLibrary.get():
            fnSmall = small.getFileName()
            fnBase = os.path.splitext(os.path.split(fnSmall)[1])[0]
            if not fnBase in smallDict:
                smallDict[fnBase]=fnSmall

        allSmallList = []
        for grid in self.inputGridSet.get():

            gridId = grid.getObjId()
            gridDir = 'grid_{}/'.format(gridId)
            fnStruct = grid.structureFile.get()
            if hasattr(grid, "bindingSiteScore"):
                bindingSiteScore = grid.bindingSiteScore.get()
            else:
                bindingSiteScore = None
            if hasattr(grid, "bindingSiteDScore"):
                bindingSiteDScore = grid.bindingSiteDScore.get()
            else:
                bindingSiteDScore = None

            smallList = []
            with open(self._getPath(gridDir + 'job_{}_pv.csv'.format(gridId))) as fhCsv:
                fnPv = self._getPath(gridDir + 'job_{}_pv.maegz'.format(gridId))
                i = 0
                for line in fhCsv.readlines():
                    if i > 1:
                        tokens = line.split(',')
                        small = SmallMolecule(smallMolFilename=smallDict[tokens[0]])
                        small.dockingScore = pwobj.Float(tokens[1])
                        small.ligandEfficiency = pwobj.Float(tokens[2])
                        small.ligandEfficiencySA = pwobj.Float(tokens[3])
                        small.ligandEfficiencyLn = pwobj.Float(tokens[4])
                        small.poseFile = pwobj.String("%d@%s"%(i, fnPv))
                        small.structFile = pwobj.String(fnStruct)
                        small.setGridId(gridId)
                        if bindingSiteScore:
                            small.bindingSiteScore = pwobj.Float(bindingSiteScore)
                        if bindingSiteDScore:
                            small.bindingSiteDScore = pwobj.Float(bindingSiteDScore)
                        smallList.append(small)
                    i += 1

            allSmallList += smallList

            mae = SchrodingerPoses(filename=fnPv)
            self._defineOutputs(**{'outputPoses_{}'.format(gridId): mae})
            self._defineSourceRelation(grid, mae)
            self._defineSourceRelation(self.inputLibrary, mae)

            if not self.mergeOutput:
                idxSorted = sortDockingResults(smallList)
                outputSet = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix=gridId)
                for idx in idxSorted:
                    small = smallList[idx]
                    outputSet.append(small)

                self._defineOutputs(**{'outputSmallMolecules_{}'.format(gridId): outputSet})
                self._defineSourceRelation(grid, outputSet)
                self._defineSourceRelation(self.inputLibrary, outputSet)
                if self.doConvertOutput:
                    print('Converting output to {}: ', 'outputSmallMolecules_{}'.format(
                        self.getEnumText('convertType'), gridId))
                    self.convertOutput(outputSet, nameDir='outputSmallMolecules_{}'.format(gridId))

        if self.mergeOutput:
            idxSorted = sortDockingResults(allSmallList)
            outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
            for idx in idxSorted:
                small = allSmallList[idx]
                outputSet.append(small)

            self._defineOutputs(outputSmallMolecules=outputSet)
            self._defineSourceRelation(self.inputGridSet, outputSet)
            self._defineSourceRelation(self.inputLibrary, outputSet)
            if self.doConvertOutput:
                print('Converting output to {}: outputSmallMolecules'.format(self.getEnumText('convertType')))
                self.convertOutput(outputSet, nameDir="outputSmallMolecules")


    def _validate(self):
        errors = []
        if self.HbondAcceptHalo.get() and self.HbondDonorHalo.get():
            errors.append('Halogens cannot be simultaneously acceptors and donors of H-bonds')
        return errors

    ###########################  Utils functions #####################
    def convertOutput(self, molSet, nameDir):
        outDir = self._getExtraPath(nameDir)
        os.mkdir(outDir)

        nt = self.numberOfThreads.get()
        nMols = len(molSet) // nt
        it, curMolSet = 0, []
        threads=[]
        for mol in molSet:
            curMolSet.append(mol.clone())
            if len(curMolSet) == nMols and it < nt - 1:
                t = threading.Thread(target=self.convertOutputStep, args=(curMolSet.copy(), outDir, it,),
                                     daemon=False)
                t.start()
                threads.append(t)
                curMolSet, it = [], it+1

        if len(curMolSet) > 0:
            t = threading.Thread(target=self.convertOutputStep, args=(curMolSet.copy(), outDir, it,),
                                 daemon=False)
            t.start()
            threads.append(t)
        for t in threads:
            t.join()

    def convertOutputStep(self, curMolSet, outDir, it):
        for i, mol in enumerate(curMolSet):
            fnAux = os.path.abspath(self._getExtraPath("tmp_%d_%d.mae" % (it, i)))
            n, fnRaw = mol.poseFile.get().split('@')
            fnOut = os.path.join(outDir, '{}.{}'.format(
                mol.getUniqueName(), self.getEnumText('convertType')))

            if not os.path.exists(fnOut):
                args = "-n %s %s -o %s" % (n, os.path.abspath(fnRaw), fnAux)
                subprocess.call([maeSubsetProg, *args.split()])

                args = fnAux + ' ' + os.path.abspath(fnOut)
                subprocess.call([structConvertProg, *args.split()])
                os.remove(fnAux)

    def convert2Mol2(self, fnSmall, it):
        fnBase = os.path.splitext(os.path.split(fnSmall)[1])[0]
        args = fnSmall
        fnSmall = self._getTmpPath('ligand{}.mol2'.format(it))
        args += " %s" % fnSmall
        subprocess.call([structConvertProg, *args.split()])
        putMol2Title(fnSmall, fnBase)
        return fnSmall

    def getAllLigandsFile(self, suffix=''):
        return 'allMoleculesFile{}.mol2'.format(suffix)
