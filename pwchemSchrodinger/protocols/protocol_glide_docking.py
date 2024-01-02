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
import os, shutil, threading, subprocess, time, glob

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam, FloatParam, IntParam, STEPS_PARALLEL
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.utils.path import makePath

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import relabelAtomsMol2, calculate_centerMass, getBaseName, performBatchThreading
from pwchem import Plugin as pwchemPlugin

from .. import Plugin as schrodinger_plugin
from ..protocols.protocol_preparation_grid import ProtSchrodingerGrid
from ..utils.utils import putMolFileTitle, convertMAEMolSet

glideProg = schrodinger_plugin.getHome('glide')
progLigPrep = schrodinger_plugin.getHome('ligprep')
structConvertProg = schrodinger_plugin.getHome('utilities/structconvert')
structCatProg = schrodinger_plugin.getHome('utilities/structcat')
propListerProg = schrodinger_plugin.getHome('utilities/proplister')
maeSubsetProg = schrodinger_plugin.getHome('utilities/maesubset')

class ProtSchrodingerGlideDocking(ProtSchrodingerGrid):
    """Calls glide to perform a docking of a set of compounds in a structural region defined by a grid.
       It is assumed that the input library of ligands is already prepared.

       The dockinsScore is the Glide Docking Score and it is measured in kcal/mol"""
    _label = 'docking (glide)'
    _program = ""

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineGlideReceptorParams(self, form):
        group = form.addGroup('Receptor')
        group.addParam('fromPockets', EnumParam, label='Dock on : ', default=1,
                       choices=['Whole protein', 'SetOfStructROIs', 'Schrodinger Grids'],
                       help='Whether to dock on a whole protein surface or on specific regions defines as StructROIs or'
                            ' Schrodinger Grids')

        group.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                       condition='fromPockets==0', label='Input structure:',
                       help='Input atomic structure to perform the docking on. We suggest the use of Structural ROIs'
                            'to speed up the docking process instead of docking on the whole protein')
        group.addParam('radius', FloatParam, label='Grid radius for whole protein: ',
                       condition='fromPockets == 0', allowsNull=False,
                       help='Radius of the Schrodinger grid for the whole protein')

        group.addParam('inputStructROIs', PointerParam, pointerClass="SetOfStructROIs",
                       condition='fromPockets==1', label='Input structural ROIs:',
                       help='Input structural ROIs defining the space where the docking will be performed')
        group.addParam('inputGridSet', PointerParam, pointerClass="SetOfSchrodingerGrids",
                       condition='fromPockets==2', label='Input grids:',
                       help='Input grids defining the space where the docking will be performed')
        return form

    def _defineGlideParams(self, form, condition='True'):
        group = form.addGroup('Docking', condition=condition)
        group.addParam('inputLibrary', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Input small molecules:', help='Input small molecules to be docked with Glide')
        group.addParam('convertOutputParam', BooleanParam, default=False, label='Convert output to mol2:',
                       help='Whether to convert the output poses from the mae to mol2 format')

        group = form.addGroup('Define grids', condition='fromPockets!=2')
        group.addParam('innerAction', EnumParam, default=1, label='Determine inner box: ',
                       choices=['Manually', 'PocketDiameter'], display=EnumParam.DISPLAY_HLIST,
                       help='How to set the inner box.'
                            'Manually: you will manually set the same x,y,z for every ROI'
                            'PocketDiameter: the diameter * n of each ROI will be used. You can set n')

        line = group.addLine('Inner box (Angstroms)', condition='innerAction==0',
                             help='The docked ligand mass center must be inside the inner box radius')
        line.addParam('innerX', IntParam, default=10, label='X')
        line.addParam('innerY', IntParam, default=10, label='Y')
        line.addParam('innerZ', IntParam, default=10, label='Z')

        group.addParam('diameterNin', FloatParam, default=0.8, condition='innerAction==1',
                       label='Size of inner box vs diameter: ',
                       help='The diameter * n of each ROI will be used as inner box side')

        group.addParam('outerAction', EnumParam, default=1, label='Determine outer box: ',
                       choices=['Manually', 'PocketDiameter'], display=EnumParam.DISPLAY_HLIST,
                       help='How to set the outer box.'
                            'Manually: you will manually set the same x,y,z for every structural ROI'
                            'PocketDiameter: the diameter * n of each pocket will be used. You can set n')

        line = group.addLine('Outer box (Angstroms)', condition='outerAction==0',
                             help='The docked ligand atoms must be inside the outer box radius.')
        line.addParam('outerX', IntParam, default=30, label='X')
        line.addParam('outerY', IntParam, default=30, label='Y')
        line.addParam('outerZ', IntParam, default=30, label='Z')

        group.addParam('diameterNout', FloatParam, default=1.2, condition='outerAction==1',
                       label='Size of outer box vs diameter: ',
                       help='The diameter * n of each structural ROI will be used as outer box side')

        group = form.addGroup('Docking')
        group.addParam('posesPerLig', IntParam, default=5, label='No. Poses to report per ligand: ',
                       help='Maximum number of final poses to report per ligand')
        group.addParam('dockingMethod', EnumParam, default=0, label='Docking method',
                       choices=['Flexible dock (confgen)', 'Rigid dock (rigid)', 'Refine (do not dock, mininplace)',
                                'Score in place (do not dock, inplace)'],
                       help='Glide method to use for docking')
        group.addParam('dockingPrecision', EnumParam, default=0, choices=['Low (HTVS)', 'Medium (SP)', 'High (XP)'],
                       label='Docking precision',
                       help='You may use a low to high strategy. HTVS takes about 2 s/ligand, SP about 10s, and XP '
                            'about 10 min.')
        group.addParam('maxkeep', IntParam, default=5000, expertLevel=LEVEL_ADVANCED,
                       label='No. Poses to keep per ligand (dock):',
                       help='Number of poses per ligand to keep in initial phase of docking.')
        group.addParam('scoreCutoff', FloatParam, default=100.0, expertLevel=LEVEL_ADVANCED,
                       label='Score cutoff:',
                       help='Scoring window for keeping initial poses.')
        group.addParam('maxref', IntParam, default=-1, expertLevel=LEVEL_ADVANCED,
                       label='No. Poses to keep per ligand (em):',
                       help='Number of poses to keep per ligand for energy minimization. '
                            'If set to -1, the default value is 400, except for XP precision that is 800.')
        group.addParam('sampleNinversions', BooleanParam, default=True, condition='dockingMethod==0',
                       label='Sample pyramid nitrogen inversions:')
        group.addParam('sampleRings', BooleanParam, default=True, condition='dockingMethod==0',
                       label='Sample rings:')
        group.addParam('epikPenalties', BooleanParam, default=False, label='Epik penalties:',
                       help='Apply penalties for ionization or tautomeric states calculated by Epik')
        group.addParam('skipMetalEpik', BooleanParam, default=True, label='Skip Epik metal only:',
                       help='Skip Epik-generated states of ligands that are designed for binding to metals. '
                            'This option is useful if the receptor has a metal but the ligand does not bind to it. '
                            'These states are skipped by default if the receptor does not have a metal.')
        group.addParam('expandedSampling', BooleanParam, default=False, label='Expanded sampling:',
                       help='Expand the sampling by bypassing the elimination of poses in the rough scoring stage. '
                            'Useful for fragment docking.')
        return form

    def _defineParams(self, form):
        form.addSection(label='Input')
        form = self._defineGlideReceptorParams(form)
        form = self._defineGridParams(form, notManualCondition='fromPockets!=2')

        group = form.addGroup('Ligands')
        group.addParam('inputLibrary', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Input small molecules:', help='Input small molecules to be docked with Glide')

        group.addParam('convertOutput2Mol2', BooleanParam, label='Convert output to mol2: ', default=False,
                       help='Whether to convert to output molecules to mol2 files instead of keeping the mae files')

        form.addSection(label='Docking')
        form = self._defineGlideParams(form)

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        convStep = self._insertFunctionStep('convertStep', prerequisites=[])
        cSteps = [convStep]

        inputGrids = self.inputGridSet.get()
        if self.fromPockets.get() == 0:
            dStep = self._insertFunctionStep('gridStep', self.inputAtomStruct.get(), prerequisites=cSteps)
            cSteps, inputGrids = [dStep], [self.inputAtomStruct.get()]

        elif self.fromPockets.get() == 1:
            gridSteps = []
            for pocket in self.inputStructROIs.get():
                dStep = self._insertFunctionStep('gridStep', pocket.clone(), prerequisites=cSteps)
                gridSteps.append(dStep)
            cSteps, inputGrids = gridSteps, self.inputStructROIs.get()

        dockSteps = []
        for grid in inputGrids:
            dStep = self._insertFunctionStep('dockingStep', grid.clone(), prerequisites=cSteps)
            dockSteps.append(dStep)
        self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

    def getLigSetFormats(self, ligSet):
        formats = []
        for lig in ligSet:
            ext = os.path.splitext(lig.getFileName())[1]
            formats.append(ext)
        formats = list(set(formats))
        return formats

    def convertStep(self):
        nt = self.numberOfThreads.get()
        # Convert receptor
        inFile = self.getOriginalReceptorFile()

        if '.mae' not in inFile:
            maeFile = self._getExtraPath('inputReceptor.maegz')
            prog = schrodinger_plugin.getHome('utilities/prepwizard')
            args = ' -WAIT -noprotassign -noimpref -noepik {} {}'. \
                format(os.path.abspath(inFile), os.path.abspath(maeFile))
            self.runJob(prog, args, cwd=self._getExtraPath())
        else:
            ext = os.path.splitext(inFile)[1]
            shutil.copy(inFile, self._getExtraPath('inputReceptor{}'.format(ext)))

        # Create ligand files
        ligSet = self.inputLibrary.get()


        ligFormats = self.getLigSetFormats(ligSet)
        if len(ligFormats) > 1 or ligFormats[0] not in ['.mae', '.sdf', '.mol2']:
            allLigandsFile = self.getAllLigandsFile(format='.mol2')
            performBatchThreading(self.createLigandsFile, ligSet, nt)
            with open(allLigandsFile, 'w') as f:
                for it in range(nt):
                    with open(self.getAllLigandsFile(suffix=it, format='.mol2')) as fLig:
                        f.write(fLig.read())
        else:
            allLigandsFile = self.getAllLigandsFile(format=ligFormats[0])
            with open(allLigandsFile, 'w') as f:
                for lig in ligSet:
                    with open(lig.getFileName()) as fLig:
                        f.write(fLig.read())

    def gridStep(self, pocket):
        if self.fromPockets.get() == 0:
            inAS = self.inputAtomStruct.get()
            pdbFile = inAS.getFileName()
            if not pdbFile.endswith('.pdb'):
                pdbFile = inAS.convert2PDB(self._getExtraPath('inputStructure.pdb'))

            pocket = self.radius.get()
            _, x, y, z = calculate_centerMass(pdbFile)
            pId = 1
        else:
            x, y, z = pocket.calculateMassCenter()
            pId = pocket.getObjId()

        fnGridDir = self.getGridDir(pId)
        gridName = self.getGridName(pId)
        # print('Grid name: ', gridName)
        makePath(fnGridDir)

        fnJob = os.path.abspath(os.path.join(fnGridDir, gridName)) + '.inp'
        fh = open(fnJob, 'w')
        fh.write("GRIDFILE %s.zip\n" % gridName)
        fh.write("OUTPUTDIR %s\n" % fnGridDir)
        fh.write("RECEP_FILE %s\n" % os.path.abspath(self.getInputMaeFile()))
        fh.write("REC_MAECHARGES True\n")
        fh.write("HBOND_DONOR_AROMH %s\n" % self.HbondDonorAromH.get())
        if self.HbondDonorAromH.get():
            fh.write("HBOND_DONOR_AROMH_CHARGE %f\n" % self.HbondDonorAromHCharge.get())
        fh.write("HBOND_ACCEP_HALO %s\n" % self.HbondAcceptHalo.get())
        fh.write("HBOND_DONOR_HALO %s\n" % self.HbondDonorHalo.get())
        fh.write("INNERBOX {},{},{}\n".format(*(self.getInnerBox(pocket))))
        fh.write("ACTXRANGE %d\n" % self.getOuterBox(pocket)[0])
        fh.write("ACTYRANGE %d\n" % self.getOuterBox(pocket)[1])
        fh.write("ACTZRANGE %d\n" % self.getOuterBox(pocket)[2])
        fh.write("OUTERBOX {},{},{}\n".format(*(self.getOuterBox(pocket))))
        fh.write("GRID_CENTER %s,%s,%s\n" % (x, y, z))
        fh.close()

        args = "-WAIT -LOCAL %s.inp" % (gridName)
        self.runJob(schrodinger_plugin.getHome('glide'), args, cwd=fnGridDir)

    def dockingStep(self, grid):
        gridId = grid.getObjId()
        if self.fromPockets.get() == 0:
            gridId = 1
            gridDir = self.getGridDir(gridId)
        elif self.fromPockets.get() == 1:
            gridDir = self.getGridDir(grid.getObjId())

        elif self.fromPockets.get() == 2:
            gridDir = self._getExtraPath('grid_{}/'.format(gridId))

            makePath(gridDir)
            fnGrid = os.path.join(gridDir, "grid_{}.zip".format(gridId))
            if not os.path.exists(fnGrid): # Prepared to resume
                shutil.copy(grid.getFileName(), fnGrid)

        fnIn = os.path.join(gridDir, 'job_{}.inp'.format(gridId))
        if not os.path.exists(fnIn): # Prepared to resume
            with open(fnIn, 'w') as fhIn:
                fhIn.write("GRIDFILE %s\n" % ("grid_{}.zip".format(gridId)))

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
                fhIn.write("HBOND_DONOR_AROMH %s\n"%self.HbondDonorAromH.get())
                if self.HbondDonorAromH.get():
                    fhIn.write("HBOND_DONOR_AROMH_CHARGE %f\n" % self.HbondDonorAromHCharge.get())
                fhIn.write("HBOND_ACCEP_HALO %s\n"%self.HbondAcceptHalo.get())
                fhIn.write("HBOND_DONOR_HALO %s\n"%self.HbondDonorHalo.get())

                fhIn.write("MAXKEEP %d\n"%self.maxkeep.get())
                fhIn.write("SCORING_CUTOFF %f\n"%self.scoreCutoff.get())
                if self.maxref.get()>0:
                    maxRefValue = self.maxref.get()
                else:
                    if self.dockingPrecision.get()==2:
                        maxRefValue = 800
                    else:
                        maxRefValue = 400
                fhIn.write("MAXREF %d\n" % maxRefValue)
                fhIn.write("POSES_PER_LIG %d\n"%self.posesPerLig.get())

                fhIn.write("LIGANDFILE {}\n".format(os.path.abspath(self.getAllLigandsFile())))

        args = "-WAIT -RESTART -LOCAL job_{}.inp".format(gridId)
        self.runJob(glideProg, args, cwd=gridDir)

        if os.path.exists(os.path.join(gridDir, "job_{}_pv.maegz".format(gridId))):
            self.runJob(propListerProg,
                        '-p "title" -p "docking score" -p "glide ligand efficiency" -p "glide ligand efficiency sa" '
                        '-p "glide ligand efficiency ln" -c -o %s %s'%\
                        ("job_{}_pv.csv".format(gridId), "job_{}_pv.maegz".format(gridId)),
                        cwd=gridDir)
        else:
            print('Failed to find ligands for grid {}'.format(gridId))

    def performOutputParsing(self, gridDirs, molLists, it, smallDict, allMaeFile):
        allSmallList = []
        for gridDir in gridDirs:
            gridId = gridDir.split('_')[1]
            gridDir = 'grid_{}/'.format(gridId)

            smallList = []
            fnPv = self._getExtraPath(gridDir + 'job_{}_pv.maegz'.format(gridId))
            with open(self._getExtraPath(gridDir + 'job_{}_pv.csv'.format(gridId))) as fhCsv:
                fhCsv.readline()
                fhCsv.readline()
                for i, line in enumerate(fhCsv.readlines()):
                    tokens = line.split(',')
                    fnBase = getBaseName(tokens[0])

                    small = SmallMolecule()
                    small.copy(smallDict[fnBase], copyId=False)
                    small._energy = pwobj.Float(tokens[1])
                    small.ligandEfficiency = pwobj.Float(tokens[2])
                    small.ligandEfficiencySA = pwobj.Float(tokens[3])
                    small.ligandEfficiencyLn = pwobj.Float(tokens[4])
                    small.setMolClass('Schrodinger')
                    small.setDockId(self.getObjId())
                    small.setGridId(gridId)

                    recFile, posFile = self.divideMaeComplex(fnPv, posIdx=i+1)
                    small.setProteinFile(recFile)
                    small.setPoseFile(posFile)
                    small.setPoseId(i + 1)

                    smallList.append(small)

            allSmallList += smallList

        if self.convertOutput2Mol2.get():
            nt = self.numberOfThreads.get()
            # print('Converting output to mol2: outputSmallMolecules')
            outDir = self._getExtraPath('outputSmallMolecules')
            os.mkdir(outDir)
            allSmallList = performBatchThreading(self.convertOutputStep, allSmallList, nt, outDir=outDir,
                                                 cloneItem=True)
        molLists[it] = allSmallList

    def createOutputStep(self):
        nt = self.numberOfThreads.get()
        smallDict = {}
        for small in self.inputLibrary.get():
            fnSmall = small.getFileName()
            fnBase = getBaseName(fnSmall)
            if fnBase not in smallDict:
                smallDict[fnBase] = small.clone()
            molName = small.getMolName()
            if molName not in smallDict:
                smallDict[molName] = small.clone()

        allMaeFile = self.mergeMAEfiles()
        fnStruct = self.getInputMaeFile()

        gridDirs = self.getGridDirs(complete=True)
        allSmallList = performBatchThreading(self.performOutputParsing, gridDirs, nt, cloneItem=False,
                                             smallDict=smallDict, allMaeFile=allMaeFile)

        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
        for small in allSmallList:
            outputSet.append(small)

        if self.convertOutputParam:
            outDir = os.path.abspath(self._getExtraPath('outputSmallMolecules'))
            os.mkdir(outDir)
            convertMAEMolSet(outputSet, outDir, nt)

        outputSet.setDocked(True)
        outputSet.proteinFile.set(self.getOriginalReceptorFile())
        outputSet.structFile = pwobj.String(fnStruct)
        if len(outputSet) > 0:
            self._defineOutputs(outputSmallMolecules=outputSet)
        else:
            print('No output docking files were generated or no poses were found')

    def divideMaeComplex(self, maeFile, posIdx=1, outDir=None):
      if not outDir:
        outDir = os.path.dirname(maeFile)
      molFile, recFile = os.path.join(outDir, getBaseName(maeFile) + f'_lig_{posIdx}.maegz'), \
                         os.path.join(outDir, getBaseName(maeFile) + '_rec.maegz')
      args = f' -n 1 {os.path.abspath(maeFile)} -o {os.path.abspath(recFile)}'
      subprocess.check_call(maeSubsetProg + args, cwd=outDir, shell=True)
      args = f' -n {posIdx+1} {os.path.abspath(maeFile)} -o {os.path.abspath(molFile)}'
      subprocess.check_call(maeSubsetProg + args, cwd=outDir, shell=True)
      return recFile, molFile

    def _validate(self):
        errors = []
        if self.HbondAcceptHalo.get() and self.HbondDonorHalo.get():
            errors.append('Halogens cannot be simultaneously acceptors and donors of H-bonds')
        return errors

    ###########################  Utils functions #####################
    def convertOutputStep(self, curMolSet, molLists, it, outDir):
        for i, mol in enumerate(curMolSet):
            molName = mol.getUniqueName()
            fnOut = os.path.join(outDir, f'{molName}.mol2')

            if not os.path.exists(fnOut):
                posFile = os.path.abspath(mol.getPoseFile())
                args = f'{posFile} {os.path.abspath(fnOut)}'
                subprocess.call([structConvertProg, *args.split()])
                fnOut = relabelAtomsMol2(fnOut)

            mol.setPoseFile(fnOut)
            molLists[it].append(mol)


    def convert2mol2(self, fnSmall, it):
        baseName = getBaseName(fnSmall)
        outFile = os.path.abspath(self._getTmpPath('{}.mol2'.format(baseName)))
        if fnSmall.endswith('.pdbqt'):
            #Manage files from autodock: 1) Convert to readable by schro (SDF). 2) correct preparation.
            # 3) Switch to mol2 to manage atom labels
            outDir = os.path.abspath(self._getTmpPath())
            args = ' -i "{}" -of sdf --outputDir "{}" --outputName {}_AD4'.format(os.path.abspath(fnSmall),
                                                                               os.path.abspath(outDir), baseName)
            pwchemPlugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=outDir, popen=True)
            auxFile = os.path.join(outDir, '{}_AD4.sdf'.format(baseName))
            fnSmall = auxFile.replace('_AD4.sdf', '_aux.sdf')
            args = " -i 0 -nt -s 1 -isd {} -osd {}".format(auxFile, fnSmall)
            if not os.path.exists(fnSmall):
                subprocess.check_call(progLigPrep + args, shell=True)

        args = " {} {}".format(os.path.abspath(fnSmall), outFile)
        if not os.path.exists(outFile):
            subprocess.check_call(structConvertProg + args, shell=True)
        while not os.path.exists(outFile):
            time.sleep(0.2)
        return outFile

    def getAllLigandsFile(self, suffix='', format=''):
        if format:
            ligFile = self._getTmpPath(f'allMoleculesFile{suffix}{format}')
        else:
            for file in os.listdir(self._getTmpPath()):
                if f'allMoleculesFile{suffix}' in file:
                    ligFile = self._getTmpPath(file)
                    break
        return ligFile

    def getGridDirs(self, complete=False):
        gridDirs = []
        for dir in os.listdir(self._getExtraPath()):
            if dir.startswith('grid_'):
                if not complete:
                    gridDirs.append(dir)
                else:
                    gridId = dir.split('_')[1]
                    if os.path.exists(self._getExtraPath(dir, "job_{}_pv.maegz".format(gridId))):
                        gridDirs.append(dir)
                    else:
                        print('No good poses found in grid ' + gridId)
        return gridDirs

    def mergeMAEfiles(self):
        maeFiles = []
        for gridDir in self.getGridDirs(complete=True):
            gridId = gridDir.split('_')[1]
            maeFiles.append(gridDir + '/job_{}_pv.maegz'.format(gridId))

        noRecMaeFiles = [maeFiles[0]]
        for mFile in maeFiles[1:]:
            mFile = os.path.abspath(self._getExtraPath(mFile))
            tFile = os.path.abspath(self._getTmpPath(getBaseName(mFile) + '.maegz'))
            args = f' -n 2: {mFile} -o {tFile}'
            self.runJob(maeSubsetProg, args, cwd=self._getTmpPath())
            noRecMaeFiles.append(tFile)

        outName = 'allMolecules.maegz'
        command = 'zcat {} | gzip -c > {}'.format(' '.join(noRecMaeFiles), outName)
        self.runJob('', command, cwd=self._getExtraPath())
        return self._getExtraPath(outName)


    def createLigandsFile(self, ligSet, molLists, it):
        curAllLigandsFile = self.getAllLigandsFile(suffix=it, format='.mol2')
        with open(curAllLigandsFile, 'w') as fh:
            for small in ligSet:
                fnSmall = small.getFileName()
                if not fnSmall.endswith('.mol2'):
                    fnSmall = self.convert2mol2(fnSmall, it)
                putMolFileTitle(fnSmall, ext='mol2')
                with open(fnSmall) as fhLigand:
                    fh.write(fhLigand.read())

    def getGridDir(self, pocketId):
        return self._getExtraPath(self.getGridName(pocketId))

    def getGridName(self, pocketId):
        return 'grid_{}'.format(pocketId)

    def getInputMaeFile(self):
        maeFile = glob.glob(self._getExtraPath('inputReceptor.mae*'))[0]
        return maeFile

    def getInnerBox(self, pocket):
        if self.fromPockets.get() == 0:
            return 3*[int(pocket * 2 * self.diameterNin.get())]
        else:
            if self.innerAction.get() == 0:
                return self.innerX.get(), self.innerY.get(), self.innerZ.get()
            else:
                diam = int(pocket.getDiameter() * self.diameterNin.get())
                return [diam, diam, diam]

    def getOuterBox(self, pocket):
        if self.fromPockets.get() == 0:
            return 3*[int(pocket * 2 * self.diameterNout.get())]
        else:
            if self.outerAction.get() == 0:
                return self.outerX.get(), self.outerY.get(), self.outerZ.get()
            else:
                diam = int(pocket.getDiameter() * self.diameterNout.get())
                return [diam, diam, diam]

    def getOriginalReceptorFile(self):
        if self.fromPockets.get() == 0:
            inAs = self.inputAtomStruct.get()
            if hasattr(inAs, '_maeFile'):
                return getattr(inAs, '_maeFile').get()
            else:
                return inAs.getFileName()

        elif self.fromPockets.get() == 1:
            inFile = self.inputStructROIs.get().getMAEFile()
            if not inFile:
                inFile = self.inputStructROIs.get().getProteinFile()

        elif self.fromPockets.get() == 2:
            inFile = self.inputGridSet.get().getProteinFile()
        return inFile



