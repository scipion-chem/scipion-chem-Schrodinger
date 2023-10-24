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
import pyworkflow.wizard as pwizard
from pwchem.wizards import AddElementSummaryWizard, DeleteElementWizard, WatchElementWizard

from ..protocols import ProtSchrodingerDesmondMD, ProtSchrodingerIFD
from ..constants import DESMOND_NPT_MD

AddElementSummaryWizard().addTarget(protocol=ProtSchrodingerDesmondMD,
                                    targets=['insertStep'],
                                    inputs=['insertStep'],
                                    outputs=['workFlowSteps', 'summarySteps'])
DeleteElementWizard().addTarget(protocol=ProtSchrodingerDesmondMD,
                                targets=['deleteStep'],
                                inputs=['deleteStep'],
                                outputs=['workFlowSteps', 'summarySteps'])
WatchElementWizard().addTarget(protocol=ProtSchrodingerDesmondMD,
                               targets=['watchStep'],
                               inputs=['watchStep'],
                               outputs=['workFlowSteps', 'summarySteps'])

AddElementSummaryWizard().addTarget(protocol=ProtSchrodingerIFD,
                                    targets=['insertStep'],
                                    inputs=['insertStep'],
                                    outputs=['workFlowSteps', 'summarySteps'])
DeleteElementWizard().addTarget(protocol=ProtSchrodingerIFD,
                                targets=['deleteStep'],
                                inputs=['deleteStep'],
                                outputs=['workFlowSteps', 'summarySteps'])
WatchElementWizard().addTarget(protocol=ProtSchrodingerIFD,
                               targets=['watchStep'],
                               inputs=['watchStep'],
                               outputs=['workFlowSteps', 'summarySteps'])

class AddDefaultStepsWizard(pwizard.Wizard):
    """Delete the step of the workflow defined by the index"""
    _targets = [(ProtSchrodingerDesmondMD, ['defSteps'])]

    def show(self, form, *params):
        protocol = form.protocol
        if protocol.defSteps.get() == protocol.DESMOND_NPT:
            form.setVar('workFlowSteps', DESMOND_NPT_MD)
            newSum = protocol.createSummary()
            form.setVar('summarySteps', newSum)
        elif protocol.defSteps.get() == protocol.NONE:
            form.setVar('workFlowSteps', '')
            newSum = protocol.createSummary()
            form.setVar('summarySteps', newSum)

