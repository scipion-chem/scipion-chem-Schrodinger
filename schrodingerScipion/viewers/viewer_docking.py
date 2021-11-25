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

from pyworkflow.protocol.params import EnumParam
from pyworkflow.viewer import DESKTOP_TKINTER
from pwchem.viewers import DockingViewer
from schrodingerScipion.protocols.protocol_glide_docking import ProtSchrodingerGlideDocking

from schrodingerScipion import Plugin
from subprocess import Popen
import os

SINGLE, MOLECULE, POCKET, MAESTRO = 'single', 'molecule', 'pocket', 'maestro'

class ProtGlideDockingViewer(DockingViewer):
    """ Visualize the output of protocol autodock """
    _label = 'Viewer glide schrodinger docking'
    _targets = [ProtSchrodingerGlideDocking]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **args):
        super().__init__(**args)
        self.maestroLabels, self.maestroLigandsDic = self.getChoices(type=MAESTRO)

    def _defineParams(self, form):
        super()._defineParams(form)
        form.addSection(label='Visualize with Maestro')
        form.addParam('displayMaestro', EnumParam,
                      choices=self.getChoices(type=MAESTRO)[0], default=0,
                      label='Display docking on Maestro: ',
                      help='Display the output docking in the Maestro software')

    def _getVisualizeDict(self):
        dispDic = super()._getVisualizeDict()
        dispDic['displayMaestro'] = self._viewMaestro
        return dispDic

    def _viewMaestro(self, e=None):
        maestroLabel = self.getEnumText('displayMaestro')
        maeFile = self.maestroLigandsDic[maestroLabel]

        cmd = [Plugin.getHome('maestro'), os.path.abspath(maeFile)]
        cwd = '/'.join(maeFile.split('/')[:-1])
        Popen(cmd, cwd=cwd, env=Plugin.getEnviron())

    def getChoices(self, type=POCKET, pymol=True):
        outputLigandsDic = {}
        for oAttr in self.protocol.iterOutputAttributes():
          if 'outputSmallMolecules' in oAttr[0] and type != MAESTRO:
            if type == POCKET:
              oLabel = oAttr[0]
            molSet = getattr(self.protocol, oAttr[0])
            for mol in molSet:
              curMol = mol.clone()
              if type == SINGLE:
                oLabel = curMol.getUniqueName()
              elif type == MOLECULE:
                oLabel = curMol.getMolBase()
              if not oLabel in outputLigandsDic:
                outputLigandsDic[oLabel] = [curMol]
              else:
                outputLigandsDic[oLabel] += [curMol]

          elif 'outputPoses' in oAttr[0] and type == MAESTRO:
              oLabel = oAttr[0]
              outputLigandsDic[oLabel] = getattr(self.protocol, oAttr[0]).getFileName()

        outputLabels = list(outputLigandsDic.keys())
        outputLabels.sort()
        if type == POCKET and pymol and len(outputLabels) > 1:
          outputLabels = ['All'] + outputLabels
        return outputLabels, outputLigandsDic

