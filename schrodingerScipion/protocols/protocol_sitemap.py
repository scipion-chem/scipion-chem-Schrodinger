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
import os
import sys

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, IntParam
from pyworkflow.utils.path import createLink
from pwem.protocols import EMProtocol
from schrodingerScipion import Plugin
from schrodingerScipion.objects import SchrodingerBindingSites

class ProtSchrodingerSiteMap(EMProtocol):
    """Calls sitemap to predict possible binding sites"""
    _label = 'binding site prediction (sitemap)'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="SchrodingerAtomStruct",
                       label='Atomic Structure:', allowsNull=False)
        form.addParam('maxsites', IntParam, expertLevel=LEVEL_ADVANCED, default=5,
                       label='Number of predicted sites:')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('sitemapStep')
        self._insertFunctionStep('createOutput')

    def sitemapStep(self):
        prog=Plugin.getHome('sitemap')

        fnIn = self._getExtraPath("atomStructIn") + self.inputStructure.get().getExtension()
        createLink(self.inputStructure.get().getFileName(), fnIn)
        fnIn = os.path.join('extra', os.path.split(fnIn)[1])

        args='-WAIT -prot %s -j job -keeplogs -keepeval'%fnIn
        args+=" -maxsites %d"%self.maxsites.get()

        self.runJob(prog,args,cwd=self._getPath())

    def createOutput(self):
        fnBinding = self._getPath("job_out.maegz")
        if os.path.exists(fnBinding):
            bindingFile=SchrodingerBindingSites(filename=fnBinding)
            bindingFile.setStructure(self.inputStructure)
            self._defineOutputs(outputGrid=bindingFile)
            self._defineSourceRelation(self.inputStructure, bindingFile)

    def _citations(self):
        return ['Halgren2009']
