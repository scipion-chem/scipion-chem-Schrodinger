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
import os, shutil
from subprocess import CalledProcessError

from pwem.objects import AtomStruct
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, IntParam, StringParam
from pyworkflow.utils.path import createLink
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from schrodingerScipion import Plugin
from schrodingerScipion.objects import SchrodingerBindingSites, SchrodingerAtomStruct
from pwchem.objects import BindingSite, SetOfBindingSites, SetOfPockets, ProteinPocket
from pwchem.constants import *
from pwchem.utils import writePDBLine, writeSurfPML, splitPDBLine

class ProtSchrodingerSiteMap(EMProtocol):
    """Calls sitemap to predict possible binding sites"""
    _label = 'binding site prediction (sitemap)'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="SchrodingerAtomStruct",
                       label='Atomic Structure:', allowsNull=False)
        form.addParam('jobName', StringParam, label='Job Name:', default='', expertLevel=LEVEL_ADVANCED)
        form.addParam('maxsites', IntParam, expertLevel=LEVEL_ADVANCED, default=5,
                       label='Number of predicted sites:')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('sitemapStep')
        self._insertFunctionStep('createOutput')

    def sitemapStep(self):
        prog=Plugin.getHome('sitemap')

        fnIn = self._getExtraPath("atomStructIn") + self.inputStructure.get().getExtension()
        createLink(self.getInputFileName(), fnIn)
        fnIn = os.path.join('extra', os.path.split(fnIn)[1])

        args='-WAIT -prot %s -j %s -keepvolpts' % (fnIn, self.getJobName())
        args+=" -maxsites %d"%self.maxsites.get()

        self.runJob(prog, args, cwd=self._getPath())

    def createOutput(self):
        fnBinding = self.getMaestroOutput()
        fnStructure = self.getInputFileName()
        fnLog = self.getOutputLogFile()
        if os.path.exists(fnBinding):
            proteinFile, pocketFiles = self.createOutputPDBFile()
            outPockets = SetOfPockets(filename=self._getPath('pockets.sqlite'))
            for oFile in pocketFiles:
              pock = ProteinPocket(os.path.abspath(oFile), os.path.abspath(proteinFile), os.path.abspath(fnLog),
                                   pClass='SiteMap')
              pock._maeFile = pwobj.String(os.path.abspath(fnStructure))
              outPockets.append(pock)

            pdbOutFile = outPockets.buildPDBhetatmFile()
            self._defineOutputs(outputPockets=outPockets)

    def _citations(self):
        return


########################## UTILS FUNCTIONS
    def getJobName(self):
      if self.jobName.get() != '':
        return self.jobName.get()
      else:
        return self.getInputFileName().split('/')[-1].split('.')[0]

    def createOutputPDBFile(self):
      maeFile = os.path.abspath(self.getMaestroOutput())
      pdbOutFile = self.getJobName()+'_out.pdb'
      pdbFiles = self.maestro2pdb(maeFile, pdbOutFile, outDir=self._getExtraPath())

      # Creates a pdb with the HETATM corresponding to pocket points
      pdbFiles, proteinFile = self.mergePDBFiles(pdbFiles, pdbOutFile)
      return proteinFile, pdbFiles

    def getMaestroOutput(self):
        return self._getPath("{}_out.maegz".format(self.getJobName()))

    def getInputFileName(self):
        return self.inputStructure.get().getFileName()

    def getOutputLogFile(self):
        return self._getPath('{}.log'.format(self.getJobName()))

    def maestro2pdb(self, maeIn, pdbOut, outDir):
      '''Convert a maestro file (.mae) to a pdb file(s)
      maeIn: input maestro file (if contains several models, there will be several outputs
      pdbOut: name of the output (with or without .pdb)'''
      pdbName, pdbOut = self.getPDBName(pdbOut)

      prog = Plugin.getHome('utilities/structconvert')
      args = '{} {}'.format(maeIn, pdbOut)
      try:
          self.runJob(prog, args, cwd=self._getExtraPath())
      except CalledProcessError as exception:
          # ask to Schrodinger why it returns a code 2 if it worked properly
          if exception.returncode != 2:
              raise exception

      pdbFiles = self.searchOutPDBFiles(pdbOut, outDir)
      if len(pdbFiles) > 1:
        return pdbFiles
      else:
        return pdbFiles[0]

    def formatPocketStrLine(self, line, numId):
      line = splitPDBLine(line)
      replacements = ['HETATM', line[1], 'APOL', 'STP', 'C', numId, *line[5:-1], '', 'Ve']
      pdbLine = writePDBLine(replacements)
      return pdbLine

    def mergePDBFiles(self, pdbFiles, pdbOutFile):
        atomLines, hetatmLines = '', ''
        idsDic = {}
        for pFile in pdbFiles:
            fileId = pFile.split('-')[1].split('.')[0]
            with open(pFile) as fpdb:
                for line in fpdb:
                    if line.startswith('TITLE'):
                      numId = line.split('_site_')[1].strip()
                      idsDic[fileId] = numId
                    elif line.startswith('ATOM'):
                      atomLines += line
                    elif line.startswith('HETATM'):
                      newLine = self.formatPocketStrLine(line, numId)
                      hetatmLines += newLine

        with open(self._getExtraPath(pdbOutFile), 'w') as f:
          f.write(atomLines)
          f.write(hetatmLines)
          f.write('\nTER')

        pdbFiles, proteinFile = self.renamePDBFiles(pdbFiles, idsDic)
        return pdbFiles, proteinFile

    def renamePDBFiles(self, pdbFiles, idsDic):
        tmpFiles = []
        for pFile in pdbFiles:
            fileId = pFile.split('-')[1].split('.')[0]
            if fileId in idsDic:
                tmpFile = pFile.replace('-{}.pdb'.format(fileId), '-{}tmp.pdb'.format(idsDic[fileId]))
                shutil.move(pFile, tmpFile)
                tmpFiles.append(tmpFile)
            else:
                proteinFile = pFile.replace('out-{}.pdb'.format(fileId), 'protein.pdb')
                shutil.move(pFile, proteinFile)

        finalPDBFiles = []
        for tmpFile in tmpFiles:
            finalPDBFiles.append(tmpFile.replace('tmp.pdb', '.pdb'))
            shutil.move(tmpFile, finalPDBFiles[-1])
        return finalPDBFiles, proteinFile

    def getPDBName(self, pdbOut):
      if '.pdb' in pdbOut:
          pdbName = pdbOut.split('.pdb')[0]
      else:
          pdbName = pdbOut[:]
          pdbOut += '.pdb'
      return pdbName, pdbOut

    def searchOutPDBFiles(self, pdbOutFile, outDir):
      pdbFiles = []
      pdbName, _ = self.getPDBName(pdbOutFile)
      for outFile in os.listdir(outDir):
        if pdbName in outFile and '.pdb' in outFile:
          pdbFiles.append(os.path.join(outDir, outFile))
      return pdbFiles








