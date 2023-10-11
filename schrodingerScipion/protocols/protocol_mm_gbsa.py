# **************************************************************************
# *
# * Authors: Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *          Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import os, time, gzip, sys

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toPdb

from pwchem.utils import pdbqt2other, convertToSdf, performBatchThreading

from ..utils import getNumberOfStructures
from .. import Plugin as schrodingerPlugin

structConvertProg = schrodingerPlugin.getHome('utilities/structconvert')
progLigPrep = schrodingerPlugin.getHome('ligprep')
progPrepWizard = schrodingerPlugin.getHome('utilities/prepwizard')
maeSubsetProg = schrodingerPlugin.getHome('utilities/maesubset')
mmgbsaProg = schrodingerPlugin.getHome('prime_mmgbsa')

class ProtSchrodingerMMGBSA(EMProtocol):
    """Optimizes the docking position and claculates the binding energy"""
    _label = 'mm-gbsa'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputSmallMolecules', params.PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Input small molecules: ',
                       help="Input docked small molecules to be optimized and scored with MM-Prime GBSA")

        group = form.addGroup('Preparation')
        group.addParam('prepareLigand', params.BooleanParam, default=True, label='Prepare ligands: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Prepare target with Schrodinger LigPrep to ensure correct structure in case it was not '
                            'originally a Schrodinger molecule. It may subtly variate the ligand atom positions. '
                            'Might be necessary if docking was not performed in Schrodinger')
        form.addParam('prepareReceptor', params.BooleanParam, default=True, label='Prepare target: ',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Prepare receptor with Schrodinger PrepWizard to ensure correct structure in case it was '
                           'not originally a Schrodinger receptor. It may subtly variate the atom positions. '
                           'Might be necessary if docking was not performed in Schrodinger')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        """ This function inserts all steps functions that will be run when running the protocol """
        convStep = self._insertFunctionStep('convertStep', prerequisites=[])
        mmStep = self._insertFunctionStep('runMMGBSAStep', prerequisites=[convStep])
        # self._insertFunctionStep('createOutputStep', prerequisites=[mmStep])
        # todo: add parse and create output

    def convertStep(self):
        nt = self.numberOfThreads.get()
        outDir = self.getInputComplexDir()
        inMols = self.inputSmallMolecules.get()
        fMol = inMols.getFirstItem()

        if '.mae' in fMol.getPoseFile():
            dockFiles = self.getMaestroDockGroups(molSet=inMols)
            for dockFile in dockFiles:
                dockFiles = self.divideInThreads(dockFile, nt, outDir)
        else:
            print('Docking files are being converted to Maestro format')
            maeFiles = performBatchThreading(self.prepareLigandFile, inMols, nt)
            recMaeFile = self.prepareReceptorFile(inMols)

            nSubMols, lastMol = self.getNMolsPerThread(len(maeFiles), nt), 0
            for i, nSubMol in enumerate(nSubMols):
                complexFile = os.path.join(outDir, os.path.basename(recMaeFile).replace('.mae', f'_{i}.mae'))
                self.runJob('zcat', '{} {} > {}'.format(recMaeFile, ' '.join(maeFiles[lastMol:lastMol+nSubMol]), complexFile),
                            cwd=self._getTmpPath())
                lastMol = lastMol+nSubMol

    def runMMGBSAStep(self):
        """ This function runs the schrodinger binary file with the given params """
        nt = self.numberOfThreads.get()
        dockDir = self.getInputComplexDir()
        dockFiles = [os.path.join(dockDir, dockFile) for dockFile in os.listdir(dockDir)]
        performBatchThreading(self.performMMGBSA, dockFiles, nt, cloneItem=False)

    # --------------------------- Parallelization functions --------------------
    def prepareLigandFile(self, mols, molLists, it):
        '''Return the corresponding MAE file for the molecule'''
        for mol in mols:
            molFile = mol.getPoseFile()
            baseName = os.path.splitext(os.path.basename(molFile))[0]
            if not molFile.endswith('.sdf'):
                # Use openbabel first to convert to sdf (usable by Schrodinger)
                molFile = convertToSdf(self, molFile)

            if self.prepareLigand.get():
                # Prepare sdf file using LigPrep
                tmpmaeFile = os.path.abspath(self._getTmpPath(baseName + '_tmp.maegz'))
                args = " -WAIT -R h -a -isd {} -omae {}".format(molFile, tmpmaeFile)
                self.runJob(progLigPrep, args, cwd=self._getTmpPath())
                maeFile = os.path.abspath(self._getTmpPath(baseName + '.maegz'))
                os.rename(tmpmaeFile, maeFile)

            else:
                # Convert sdf to mae file
                maeFile = os.path.abspath(self._getTmpPath(mol.getUniqueName() + '.maegz'))
                self.runJob(structConvertProg, '{} {}'.format(molFile, maeFile))

            molLists[it].append(maeFile)

    def performMMGBSA(self, dockFiles, molLists, it):
        #  todo: add args
        for dockFile in dockFiles:
            print('Performing prime mm-gbsa on : ', dockFile)
            args = f"-WAIT {os.path.abspath(dockFile)}"
            self.runJob(mmgbsaProg, args, cwd=self._getExtraPath())

    # --------------------------- Utils functions --------------------
    def getNMolsPerThread(self, nMols, nt):
        '''Return a list with the number of element to process per thread'''
        nSubMols = [0 for _ in range(nt)]
        for i in range(nMols):
            nSubMols[i % nt] += 1
        return nSubMols

    def divideInThreads(self, dockFile, nt, outDir):
        '''Divide the complex file into several containing the same receptor (first mol) so the following
        process of these complex files can be parallelized in nt threads'''
        nMols = getNumberOfStructures(dockFile) - 1
        if nMols < nt:
            nt = nMols

        if nt > 1:
            nSubMols = self.getNMolsPerThread(nMols, nt)
            firstMol, maeFiles = 2, []
            for i, nSubMol in enumerate(nSubMols):
                maeFile = os.path.basename(dockFile).replace('.mae', f'_{i}.mae')
                lastMol = firstMol + nSubMol
                args = ' -n "1, {}:{}" {} -o {}'.format(firstMol, lastMol, os.path.abspath(dockFile), maeFile)
                self.runJob(maeSubsetProg, args, cwd=outDir)
                firstMol = lastMol
                maeFiles.append(os.path.join(outDir, maeFile))
        else:
            maeFiles = [os.path.join(outDir, os.path.basename(dockFile))]
            os.link(dockFile, maeFiles[0])
        return maeFiles

    def getMaestroDockGroups(self, molSet):
        '''Return the complex (receptor-molecules) files in a molSet'''
        gFiles = []
        for mol in molSet:
            g = mol.getPoseFile().split('@')[1]
            if not g in gFiles:
                gFiles.append(g)
        return gFiles

    def getInputComplexDir(self):
        cDir = self._getExtraPath('inputComplex')
        if not os.path.exists(cDir):
            os.mkdir(cDir)
        return os.path.abspath(cDir)

    def getInputComplexFiles(self):
        cFiles = []
        cDir = self.getInputComplexDir()
        for file in os.listdir(cDir):
            cFiles.append(os.path.join(cDir, file))
        return cFiles

    def getPdbFile(self, molSet):
      proteinFile = molSet.getProteinFile()
      inName, inExt = os.path.splitext(os.path.basename(proteinFile))

      if inExt == '.pdb':
          return os.path.abspath(proteinFile)
      else:
        pdbFile = os.path.abspath(os.path.join(self._getTmpPath(inName + '.pdb')))
        if inExt == '.pdbqt':
            pdbqt2other(self, proteinFile, pdbFile)
        else:
            toPdb(proteinFile, pdbFile)
        return os.path.abspath(pdbFile)

    def prepareReceptorFile(self, molSet):
        '''Return the corresponding MAE file for the receptor the molecule is docked to'''
        mol = molSet.getFirstItem()
        if hasattr(mol, 'structFile'):
            targetMaeFile = molSet.structFile
        else:
            pdbFile = self.getPdbFile(molSet)
            targetName = os.path.splitext(os.path.basename(pdbFile))[0]
            targetMaeFile = os.path.abspath(self._getTmpPath(targetName + '.maegz'))
            if self.prepareReceptor.get():
                args = '-WAIT -noprotassign -noimpref -noepik '
                args += '%s %s' % (os.path.abspath(pdbFile), os.path.abspath(targetMaeFile))
                self.runJob(progPrepWizard, args, cwd=self._getTmpPath())
            else:
                self.runJob(structConvertProg, '{} {}'.format(pdbFile, targetMaeFile), cwd=self._getTmpPath())
        return targetMaeFile

