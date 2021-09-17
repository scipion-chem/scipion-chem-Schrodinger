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

VOLUME_PYMOL, VOLUME_PYMOL_SURF = 1, 2

class ViewerSitemap(pwviewer.ProtocolViewer):
  _label = 'Viewer pockets'
  _targets = [ProtSchrodingerSiteMap]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Visualization of consensus pockets')
    form.addParam('output', params.EnumParam,
                  choices=['SetOfPockets', 'AtomStructure', 'SchrodingerAtomStructure'],
                  default=0, display=params.EnumParam.DISPLAY_HLIST,
                  label='Output to display: ',
                  help='Scipion output to visualize'
                  )

    form.addParam('displayOutput', params.LabelParam,
                  condition='output!=0',
                  label='Display output: ',
                  help='*AtomStruct*: display output AtomStruct using Chimera.\n '
                       '*SchrodingerAtomStructure*: display output using Maestro.'
                  )

    form.addParam('displayPyMol', params.EnumParam,
                  choices=['PyMol (Pocket Points)', 'PyMol (Contact Surface)'],
                  default=0, condition='output==0',
                  display=params.EnumParam.DISPLAY_HLIST,
                  label='Display output: ',
                  help='*PyMol*: display Set of Pockets and pockets as points / surface'
                  )

  def _getVisualizeDict(self):
    return {
      'displayOutput': self._showOutput,
      'displayPyMol': self._showStandardPockets,
    }

  def _validate(self):
    return []

  # =========================================================================
  # ShowAtomStructs
  # =========================================================================

  def _showOutput(self, paramName=None):
    if self.output.get() == 1:
      return self._showAtomStructChimera()

    elif self.output.get() == 2:
      return self._showSchAtomStructMaestro()


  def _showStandardPockets(self, paramName=None):
    if self.displayPyMol == VOLUME_PYMOL-1:
      return self._showAtomStructPyMolPoints()

    elif self.displayPyMol == VOLUME_PYMOL_SURF-1:
      return self._showAtomStructPyMolSurf()

  #Display functions
  def _showAtomStructPyMolPoints(self):
    outPockets = getattr(self.protocol, 'outputPockets')
    pymolV = PocketPointsViewer(project=self.getProject())
    pymolV._visualize(outPockets)

  def _showAtomStructPyMolSurf(self):
    outPockets = getattr(self.protocol, 'outputPockets')
    pymolV = ContactSurfaceViewer(project=self.getProject())
    pymolV._visualize(outPockets)

  def _showAtomStructChimera(self):
    outAtomStruct = getattr(self.protocol, 'outputAtomStruct')
    ChimeraView(outAtomStruct.getFileName()).show()

  def _showSchAtomStructMaestro(self):
    outSchAtomStruct = getattr(self.protocol, 'outputSchrodingerAtomStruct')
    SchV = SchrodingerDataViewer(project=self.getProject())
    SchV._visualize(outSchAtomStruct)
