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


from pyworkflow.viewer import DESKTOP_TKINTER, ProtocolViewer
from pyworkflow.protocol.params import IntParam
import pyworkflow.utils as pwutils

from schrodingerScipion import Plugin as schrodinger_plugin
from schrodingerScipion.protocols.protocol_glide_docking import ProtSchrodingerGlideDocking

class ProtSchrodingerGlideDockingViewer(ProtocolViewer):
    """ Visualize the output of protocol Glide Docking """
    _label = 'viewer glide docking'
    _targets = [ProtSchrodingerGlideDocking]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **args):
        ProtocolViewer.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('pose', IntParam, default=1,
                      label="Pose Id", help='As shown in the output SetOfSmallMolecules')

    def _getVisualizeDict(self):
        return {'pose': self._viewResults}

    def _viewResults(self, e=None):
        views = []
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
        return views
