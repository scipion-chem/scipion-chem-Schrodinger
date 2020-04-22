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
import glob
import os

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, IntParam
from pyworkflow.utils.path import createLink
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from schrodingerScipion import Plugin
from schrodingerScipion.objects import SchrodingerBindingSites
from bioinformatics.objects import BindingSite, SetOfBindingSites

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
        def parseEvalLog(fnLog):
            fh=open(fnLog)
            state = 0
            for line in fh.readlines():
                if state==0 and line.startswith("SiteScore"):
                    state=1
                elif state==1:
                    # SiteScore size   Dscore  volume  exposure enclosure contact  phobic   philic   balance  don/acc
                    tokens = [float(x) for x in line.split()]
                    return tokens
            return None

        fnBinding = self._getPath("job_out.maegz")
        if os.path.exists(fnBinding):
            setOfBindings = SetOfBindingSites().create(path=self._getPath())
            for fn in glob.glob(self._getPath("job_site_*_eval.log")):
                score, size, dscore, volume, exposure, enclosure, contact, phobic, philic, balance, donacc = parseEvalLog(fn)
                n = fn.split("job_site_")[1].replace("_eval.log","")
                bindingSite = BindingSite(bindingSiteFilename="%s@%s"%(n,fnBinding))
                bindingSite.score     = pwobj.Float(score)
                bindingSite.size      = pwobj.Float(size)
                bindingSite.dscore    = pwobj.Float(dscore)
                bindingSite.volume    = pwobj.Float(volume)
                bindingSite.exposure  = pwobj.Float(exposure)
                bindingSite.enclosure = pwobj.Float(enclosure)
                bindingSite.contact   = pwobj.Float(contact)
                bindingSite.phobic    = pwobj.Float(phobic)
                bindingSite.phobic    = pwobj.Float(phobic)
                bindingSite.philic    = pwobj.Float(philic)
                bindingSite.balance   = pwobj.Float(balance)
                bindingSite.donacc    = pwobj.Float(donacc)
                bindingSite.setStructure(self.inputStructure)

                setOfBindings.append(bindingSite)

            self._defineOutputs(outputSetBindingSites=setOfBindings)
            self._defineSourceRelation(self.inputStructure, setOfBindings)

            mae = SchrodingerBindingSites(filename=fnBinding)
            self._defineOutputs(outputBindingSites=mae)
            self._defineSourceRelation(self.inputStructure, mae)

    def _citations(self):
        return ['Halgren2009']
