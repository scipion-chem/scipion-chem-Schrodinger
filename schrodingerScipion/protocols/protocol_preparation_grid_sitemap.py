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

from pyworkflow.protocol.params import MultiPointerParam, BooleanParam, FloatParam, IntParam
from pyworkflow.object import String, Float
from pyworkflow.utils.path import createLink, makePath

from pwem.protocols import EMProtocol
from schrodingerScipion.objects import SchrodingerGrid

from schrodingerScipion import Plugin as schrodinger_plugin

class ProtSchrodingerGridSiteMap(EMProtocol):
    """Calls glide to prepare a grid with an input that is a set of binding sites"""
    _label = 'grid definition binding site (glide)'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSites', MultiPointerParam, pointerClass="BindingSite",
                       label='Binding sites:', allowsNull=False)
        line = form.addLine('Inner box (Angstroms)')
        line.addParam('innerX', IntParam, default=10, label='X')
        line.addParam('innerY', IntParam, default=10, label='Y')
        line.addParam('innerZ', IntParam, default=10, label='Z')

        line = form.addLine('Outer box (Angstroms)')
        line.addParam('outerX', IntParam, default=30, label='X')
        line.addParam('outerY', IntParam, default=30, label='Y')
        line.addParam('outerZ', IntParam, default=30, label='Z')

        form.addParam('HbondDonorAromH', BooleanParam, default=False, label='Aromatic H as H-bond donors:',
                      help='Accept aromatic hydrogens as potential H-bond donors.')
        form.addParam('HbondDonorAromHCharge', FloatParam, default=0.0, label='Aromatic H as H-bond donors Charge:',
                      condition='HbondDonorAromH',
                      help='Partial charge cutoff for accepting aromatic hydrogens as potential H-bond donors. '
                           'The cutoff is applied to the actual (signed) charge, not the absolute value.')
        form.addParam('HbondAcceptHalo', BooleanParam, default=False, label='Halogens as H-bond acceptors:',
                      help='Accept halogens (neutral or charged, F, Cl, Br, or I) as H-bond acceptors.')
        form.addParam('HbondDonorHalo', BooleanParam, default=False, label='Halogens as H-bond donors:',
                      help='Accept the halogens (Cl, Br, I, but not F) as potential H-bond '
                           '(noncovalent interaction) donors')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('preparationStep')
        self._insertFunctionStep('createOutput')

    def preparationStep(self):
        for site in self.inputSites:
            fnSite = site.get().getFileName()
            n = fnSite.split('@')[0]

            args = schrodinger_plugin.getPluginHome('utils/schrodingerUtils.py') + " centroid %s %s" %\
                   (fnSite,self._getTmpPath("centroid.txt"))
            schrodinger_plugin.runSchrodinger(self, "python3", args)
            fh = open(self._getTmpPath("centroid.txt"))
            line = fh.readline()
            fh.close()
            x,y,z = line.split()

            fnGridDir = "grid-%s"%n
            makePath(self._getPath(fnGridDir))
            fnTarget = site.get().structureFile.get()
            fnTargetLocal = self._getPath("%s/%s.maegz"%(fnGridDir,fnGridDir))
            createLink(fnTarget, fnTargetLocal)

            fnJob = self._getPath('%s/%s.inp'%(fnGridDir,fnGridDir))
            fh = open(fnJob,'w')
            fh.write("GRIDFILE %s.zip\n"%fnGridDir)
            fh.write("OUTPUTDIR %s\n"%fnGridDir)
            fh.write("RECEP_FILE %s.maegz\n"%fnGridDir)
            fh.write("REC_MAECHARGES True\n")
            fh.write("HBOND_DONOR_AROMH %s\n"%self.HbondDonorAromH.get())
            if self.HbondDonorAromH.get():
                fh.write("HBOND_DONOR_AROMH_CHARGE %f\n" % self.HbondDonorAromHCharge.get())
            fh.write("HBOND_ACCEP_HALO %s\n"%self.HbondAcceptHalo.get())
            fh.write("HBOND_DONOR_HALO %s\n"%self.HbondDonorHalo.get())
            fh.write("INNERBOX %d,%d,%d\n"%(self.innerX.get(),self.innerY.get(),self.innerZ.get()))
            fh.write("ACTXRANGE %d\n"%self.outerX.get())
            fh.write("ACTYRANGE %d\n"%self.outerY.get())
            fh.write("ACTZRANGE %d\n"%self.outerZ.get())
            fh.write("OUTERBOX %d,%d,%d\n"%(self.outerX.get(),self.outerY.get(),self.outerZ.get()))
            fh.write("GRID_CENTER %s,%s,%s\n"%(x,y,z))
            fh.close()

            args = "-WAIT -LOCAL %s.inp"%fnGridDir
            self.runJob(schrodinger_plugin.getHome('glide'), args, cwd=self._getPath(fnGridDir))

    def createOutput(self):
        for site in self.inputSites:
            fnSite = site.get().getFileName()
            n = fnSite.split('@')[0]

            fnDir = self._getPath("grid-%s"%n)
            if os.path.exists(fnDir):
                fnBase = os.path.split(fnDir)[1]
                fnGrid = os.path.join(fnDir,'%s.zip'%fnBase)
                if os.path.exists(fnGrid):
                    gridFile=SchrodingerGrid(filename=fnGrid)
                    gridFile.structureFile=String(site.get().structureFile.get())
                    gridFile.bindingSiteScore=Float(site.get().score.get())
                    gridFile.bindingSiteDScore=Float(site.get().dscore.get())

                    n = fnDir.split('grid-')[1]
                    outputDict = {'outputGrid%s' % n: gridFile}
                    self._defineOutputs(**outputDict)
                    self._defineSourceRelation(site,gridFile)

    def _validate(self):
        errors = []
        if self.HbondAcceptHalo.get() and self.HbondDonorHalo.get():
            errors.append('Halogens cannot be simultaneously acceptors and donors of H-bonds')
        return errors
