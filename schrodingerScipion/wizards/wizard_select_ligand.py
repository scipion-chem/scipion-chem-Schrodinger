# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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

"""
This wizard will extract the chains from a atomic structure (pdb) file in
order to select it in the protocol.
Then, it will load the structure and will take all chain related
information such as name and number of residues.
"""

# Imports
from pwem.wizards.wizard import EmWizard
from schrodingerScipion.protocols.protocol_desmond_systemPrep import ProtSchrodingerDesmondSysPrep

import pyworkflow.wizard as pwizard
import pyworkflow.object as pwobj
from pyworkflow.gui.tree import ListTreeProviderString
from pyworkflow.gui import dialog

class GetLigandsWizard(pwizard.Wizard):
    """Load an atomic structure, parse chain related information as
       name, number of residues, list of aminoacids (or other residues)"""
    _targets = [(ProtSchrodingerDesmondSysPrep, ['inputLigand'])]

    def getListOfMols(self, protocol):
        molsList = []
        if hasattr(protocol, 'inputSetOfMols') and protocol.inputSetOfMols.get() is not None:
            for mol in protocol.inputSetOfMols.get():
                molsList.append(mol.getUniqueName())
        return molsList

    def show(self, form, *params):
        protocol = form.protocol
        try:
            listOfMols = self.getListOfMols(protocol)
        except Exception as e:
            print("ERROR: ", e)
            return

        finalMolsList = []
        for i in listOfMols:
            finalMolsList.append(pwobj.String(i))
        provider = ListTreeProviderString(finalMolsList)
        dlg = dialog.ListDialog(form.root, "Small Molecules", provider,
                                "Select one of molecules (Pocket_MolName-Conformer_Position)")
        form.setVar('inputLigand', dlg.values[0].get())