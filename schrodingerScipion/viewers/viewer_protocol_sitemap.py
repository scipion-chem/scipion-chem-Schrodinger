# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
from ..protocols import ProtSchrodingerSiteMap
import pyworkflow.protocol.params as params
import pyworkflow.viewer as pwviewer
from pwchem.viewers import PyMolViewer, PocketPointsViewer, ContactSurfaceViewer, VmdViewFpocket
from pwem.viewers import ChimeraView
from ..viewers import SchrodingerDataViewer

VOLUME_PYMOL, VOLUME_PYMOL_SURF = 0, 1

class ViewerSitemap(pwviewer.ProtocolViewer):
  _label = 'Viewer pockets'
  _targets = [ProtSchrodingerSiteMap]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Visualization of consensus pockets')
    form.addParam('displayPyMol', params.EnumParam,
                  choices=['PyMol (Pocket Points)', 'PyMol (Contact Surface)'],
                  default=0,
                  display=params.EnumParam.DISPLAY_HLIST,
                  label='Display output: ',
                  help='*PyMol*: display Set of Pockets and pockets as points / surface'
                  )
    form.addParam('displayBBoxes', params.BooleanParam,
                  default=False, label='Display pocket bounding boxes',
                  help='Display the bounding boxes in pymol to check the size for the localized docking')
    form.addParam('pocketRadiusN', params.FloatParam, label='Grid radius vs pocket radius: ',
                  default=1.1, condition='displayBBoxes',
                  help='The radius * n of each pocket will be used as grid radius')

  def _getVisualizeDict(self):
    return {
      'displayPyMol': self._showStandardPockets,
    }

  def _validate(self):
    return []

  # =========================================================================
  # ShowAtomStructs
  # =========================================================================

  def _showStandardPockets(self, paramName=None):
    if self.displayPyMol == VOLUME_PYMOL:
      return self._showAtomStructPyMolPoints()

    elif self.displayPyMol == VOLUME_PYMOL_SURF:
      return self._showAtomStructPyMolSurf()

  #Display functions
  def _showAtomStructPyMolPoints(self):
    bBox = self.displayBBoxes.get()
    if bBox:
      bBox = self.pocketRadiusN.get()

    outPockets = getattr(self.protocol, 'outputPockets')
    pymolV = PocketPointsViewer(project=self.getProject())
    pymolV._visualize(outPockets, bBox=bBox)

  def _showAtomStructPyMolSurf(self):
    bBox = self.displayBBoxes.get()
    if bBox:
      bBox = self.pocketRadiusN.get()

    outPockets = getattr(self.protocol, 'outputPockets')
    pymolV = ContactSurfaceViewer(project=self.getProject())
    pymolV._visualize(outPockets, bBox=bBox)
