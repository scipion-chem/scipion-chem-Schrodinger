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

from pyworkflow.protocol.params import EnumParam
from pwchem.viewers import SmallMoleculesViewer

from schrodingerScipion.protocols.protocol_glide_docking import ProtSchrodingerGlideDocking
from schrodingerScipion.viewers.viewers_data import MaestroView

SINGLE, MOLECULE, POCKET, SET = 'single', 'molecule', 'pocket', 'set'

class ProtGlideDockingViewer(SmallMoleculesViewer):
    """ Visualize the output of protocol autodock """
    _label = 'Viewer glide schrodinger docking'
    _targets = [ProtSchrodingerGlideDocking]

    def __init__(self, **args):
        super().__init__(**args)

    def _defineParams(self, form):
        super()._defineParams(form)
        form.addSection(label='Maestro view')
        form.addParam('displaymaestroPocket', EnumParam,
                       choices=self.getChoices(vType=POCKET, pymol=False)[0], default=0,
                       label='Display one ligand type: ',
                       help='Display all conformers and positions of this molecule')

    def _getVisualizeDict(self):
        visDic = super()._getVisualizeDict()
        visDic.update({'displayMaestroPocket': self._viewPocketMaestroDock})
        return visDic

    def _viewPocketMaestroDock(self, e=None):
        ligandLabel = self.getEnumText('displaymaestroPocket')
        mols = self.pocketLigandsDic[ligandLabel]

        return [MaestroView(os.path.abspath(mols[0].maeFile.get()), cwd=self.protocol._getExtraPath())]
