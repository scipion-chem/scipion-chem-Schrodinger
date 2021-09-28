# **************************************************************************
# *
# * Authors:  Carlos Oscar Sorzano (coss@cnb.csic.es)
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
from subprocess import call, Popen
from pyworkflow.viewer import DESKTOP_TKINTER, ProtocolViewer
from pyworkflow.protocol.params import IntParam, EnumParam, BooleanParam
import pyworkflow.utils as pwutils
from pwchem.viewers import PyMolViewer

from schrodingerScipion import Plugin as schrodinger_plugin
from schrodingerScipion.protocols.protocol_glide_docking import ProtSchrodingerGlideDocking

PYMOL, MAESTRO = 0, 1

class ProtSchrodingerGlideDockingViewer(ProtocolViewer):
    """ Visualize the output of protocol Glide Docking """
    _label = 'viewer glide docking'
    _targets = [ProtSchrodingerGlideDocking]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **args):
        ProtocolViewer.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Docking visualization')
        form.addParam('displaySoft', EnumParam,
                      choices=['PyMol', 'Maestro'], default=PYMOL, display=EnumParam.DISPLAY_HLIST,
                      label='Display docking output with: ',
                      help='*PyMol*: display target and ligands in pymol.\n '
                           '*Maestro*: display target and ligands in Maestro.')
        form.addParam('singleLigand', BooleanParam, condition='displaySoft=={}'.format(PYMOL),
                      default=False,
                      label='Display single ligand: ',
                      help='Display the target with a single ligand docked')
        form.addParam('displayPymol', EnumParam, condition='displaySoft=={} and not singleLigand'.format(PYMOL),
                      choices=self.getChoices(pymol=True), default=0,
                      label='Display docking on pocket result: ',
                      help='Docking results are grouped by their pocket, choose the one to visualize')
        form.addParam('displayPymolSingle', EnumParam, condition='displaySoft=={} and singleLigand'.format(PYMOL),
                      choices=self.getChoicesSingle(), default=0,
                      label='Display single ligand: ',
                      help='Display this single ligand with the target')

        form.addParam('displayMaestro', EnumParam, condition='displaySoft=={}'.format(MAESTRO),
                      choices=self.getChoices(pymol=False), default=0,
                      label='Display docking on pocket result: ',
                      help='Docking results are grouped by their pocket, choose the one to visualize')


    def getChoicesSingle(self):
        self.outputLigands = {}
        for oAttr in self.protocol.iterOutputAttributes():
            if 'outputSmallMolecules' in oAttr[0]:
                molSet = getattr(self.protocol, oAttr[0])
                for mol in molSet:
                    curMol = mol.clone()
                    pName = curMol.getUniqueName()
                    self.outputLigands[pName] = curMol

        outputLabels = list(self.outputLigands.keys())
        outputLabels.sort()
        return outputLabels

    def getChoices(self, pymol=True):
        outputLabels = []
        for oAttr in self.protocol.iterOutputAttributes():
            if pymol and 'outputSmallMolecules' in oAttr[0]:
                outputLabels.append(oAttr[0])
            elif not pymol and 'outputPoses' in oAttr[0]:
                outputLabels.append(oAttr[0])
        outputLabels.sort()
        if pymol and not self.protocol.mergeOutput:
            outputLabels = ['All'] + outputLabels
        return outputLabels


    def _getVisualizeDict(self):
        return {'displayPymol': self._viewPymol,
                'displayMaestro': self._viewMaestro,
                'displayPymolSingle': self._viewSinglePymol}

    def _viewSinglePymol(self, e=None):
        #todo: develop a viewer for single poses
        ligandLabel = self.getEnumText('displayPymolSingle')
        mol = self.outputLigands[ligandLabel].clone()
        pmlFile = self.protocol._getExtraPath(ligandLabel)+'/{}.pml'.format(ligandLabel)
        self.writePmlFile(pmlFile, self.buildPMLDockingSingleStr(mol, ligandLabel, addTarget=True))

        pymolV = PyMolViewer(project=self.getProject())
        pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

    def _viewMaestro(self, e=None):
        fnPv = getattr(self.protocol, self.getEnumText('displayMaestro')).getFileName()
        pwutils.runJob(None, schrodinger_plugin.getHome('maestro'), fnPv,
                       env=schrodinger_plugin.getEnviron())

    def _viewPymol(self, e=None):
        if self.getEnumText('displayPymol') != 'All':
            outName = self.getEnumText('displayPymol')
            pmlFile = self.protocol._getPath('{}.pml'.format(outName))
            pmlStr = self.buildPMLDockingsStr(outName)
        else:
            pmlFile = self.protocol._getPath('allOutputMols.pml')
            outName = self.getChoices()[1]
            pmlStr = self.buildPMLDockingsStr(outName)
            for outName in self.getChoices()[2:]:
                pmlStr += self.buildPMLDockingsStr(outName, addTarget=False)

        self.writePmlFile(pmlFile, pmlStr)
        pymolV = PyMolViewer(project=self.getProject())
        pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))


    def buildPMLDockingSingleStr(self, mol, molName, addTarget=True):
        pmlStr = ''
        if addTarget:
            pmlStr = 'load {}\n'.format(self.protocol.inputGridSet.get().getProteinFile())

        pdbFile = os.path.abspath(self.maePos2pdb(mol, 'single', molName))
        pmlStr += 'load {}, {}\n'.format(pdbFile, molName)
        return pmlStr

    def buildPMLDockingsStr(self, outName, addTarget=True):
        molSet = getattr(self.protocol, outName)
        pmlStr = ''
        if addTarget:
            pmlStr = 'load {}\n'.format(self.protocol.inputGridSet.get().getProteinFile())
        if not os.path.exists(self.getPDBsDir(outName)):
            molDic = {}
            for i, mol in enumerate(molSet):
                pdbFile = self.maePos2pdb(mol, i, outName)
                pFile = os.path.abspath(pdbFile)
                pName = mol.getUniqueName()
                molDic[pName] = pFile

            molKeys = list(molDic.keys())
            molKeys.sort()
            for molName in molKeys:
                pFile = molDic[molName]
                pmlStr += 'load {}, {}\ndisable {}\n'.format(pFile, molName, molName)

        else:
            molFiles = list(os.listdir(self.getPDBsDir(outName)))
            molFiles.sort()
            for pdbFile in molFiles:
                pFile = os.path.abspath(os.path.join(self.getPDBsDir(outName), pdbFile))
                pName = pdbFile.split('.')[0]
                pmlStr += 'load {}, {}\ndisable {}\n'.format(pFile, pName, pName)
        return pmlStr

    def writePmlFile(self, pmlFile, pmlStr):
        with open(pmlFile, 'w') as f:
            f.write(pmlStr)
            f.write('zoom')
        return pmlFile, pmlStr

    def maePos2pdb(self, mol, i, outName):
        progMaeSubset = schrodinger_plugin.getHome('utilities/maesubset')
        progStructConvert = schrodinger_plugin.getHome('utilities/structconvert')

        outDir = self.getPDBsDir(outName)
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        fnAux = os.path.abspath(self.protocol._getExtraPath("tmp_{}.mae".format(i)))
        n, fnRaw = mol.poseFile.get().split('@')
        fnOut = os.path.join(outDir, '{}.pdb'.format(mol.getUniqueName()))

        if not os.path.basename(fnOut) in os.listdir(outDir):
            args = "-n %s %s -o %s" % (n, os.path.abspath(fnRaw), fnAux)
            p1= Popen([progMaeSubset, *args.split()])
            p1.wait()

            args = [fnAux, os.path.abspath(fnOut)]
            p2 = Popen([progStructConvert, *args])
            p2.wait()
            os.remove(fnAux)
        return fnOut

    def getPDBsDir(self, outName):
        return self.protocol._getExtraPath(outName)

