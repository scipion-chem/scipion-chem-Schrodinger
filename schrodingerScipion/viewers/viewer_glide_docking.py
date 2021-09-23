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
from pyworkflow.protocol.params import IntParam, EnumParam
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
        form.addParam('displayPymol', EnumParam, condition='displaySoft=={}'.format(PYMOL),
                      choices=self.getChoices(pymol=True), default=0,
                      label='Display docking on pocket result: ',
                      help='Docking results are grouped by their pocket, choose the one to visualize')

        form.addParam('displayMaestro', EnumParam, condition='displaySoft=={}'.format(MAESTRO),
                      choices=self.getChoices(pymol=False), default=0,
                      label='Display docking on pocket result: ',
                      help='Docking results are grouped by their pocket, choose the one to visualize')

    def _getVisualizeDict(self):
        return {'displayPymol': self._viewPymol,
                'displayMaestro': self._viewMaestro}

    def _viewPose(self, e=None):
        #todo: develop a viewer for single poses
        i = self.pose.get()
        fnPose = None
        for small in self.protocol.outputSmallMolecules:
            if small.getObjId() == i:
                fnPose = small.poseFile.get()
                break
        if fnPose:
            i = int(fnPose.split('@')[0])
            fnPv = self.protocol.outputPoses.getFileName()
            fnPose = self.protocol._getTmpPath('posei.maegz')
            pwutils.runJob(None, schrodinger_plugin.getHome('utilities/maesubset'),
                           "-n %d:%d %s -o %s" % (i, i, fnPv, fnPose))
            fnStruct = self.protocol.inputGrid.get().structureFile.get()
            fnBoth = self.protocol._getTmpPath('pv.maegz')
            pwutils.runJob(None, schrodinger_plugin.getHome('utilities/structcat'),
                           "-imae %s -imae %s -o %s" % (fnPose,fnStruct,fnBoth))

            pwutils.runJob(None, schrodinger_plugin.getHome('maestro'), fnBoth,
                           env=schrodinger_plugin.getEnviron())

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
            pmlStr = self.buildPMLDockingsStr(outName, prefix='grid-{}_'.format(outName.split('_')[-1]))
            for outName in self.getChoices()[2:]:
                pmlStr += self.buildPMLDockingsStr(outName, addTarget=False,
                                                   prefix='grid-{}_'.format(outName.split('_')[-1]))

        self.writePmlFile(pmlFile, pmlStr)
        pymolV = PyMolViewer(project=self.getProject())
        pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

    def getChoices(self, pymol=True):
        outputLabels = []
        for oAttr in self.protocol.iterOutputAttributes():
            if pymol and 'outputSmallMolecules' in oAttr[0]:
                outputLabels.append(oAttr[0])
            elif not pymol and 'outputPoses' in oAttr[0]:
                outputLabels.append(oAttr[0])
        outputLabels.sort()
        if pymol:
            outputLabels = ['All'] + outputLabels
        return outputLabels

    def buildPMLDockingsStr(self, outName, addTarget=True, prefix=''):
        molSet = getattr(self.protocol, outName)
        if addTarget:
            pmlStr = 'load {}\n'.format(self.protocol.inputGridSet.get().getProteinFile())
        else:
            pmlStr = ''
        if not os.path.exists(self.getPDBsDir(outName)):
            for i, mol in enumerate(molSet):
                pdbFile = self.maePos2pdb(mol, i, outName)
                pFile = os.path.abspath(pdbFile), os.path.abspath(pdbFile)
                pName = prefix + pdbFile.split('.')[0]
                pmlStr += 'load {}, {}\ndisable {}\n'.format(pFile, pName, pName)

        else:
            for pdbFile in os.listdir(self.getPDBsDir(outName)):
                pFile = os.path.abspath(os.path.join(self.getPDBsDir(outName), pdbFile))
                pName = prefix + pdbFile.split('.')[0]
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
        fnAux = os.path.abspath(self.protocol._getExtraPath("tmp_%d.mae" % i))
        n, fnRaw = mol.poseFile.get().split('@')
        fnOut = os.path.join(outDir, os.path.basename(mol.getFileName()).split('.')[0] + '_{}.pdb'.format(n))

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

