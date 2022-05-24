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

from pyworkflow.protocol.params import MultiPointerParam, FloatParam, IntParam, PointerParam
from pyworkflow.object import String, Float
from pyworkflow.utils.path import createLink, makePath
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pwem.protocols import EMProtocol

from pwchem.objects import BindingSite, SetOfBindingSites
from pwchem import Plugin as pwchem_plugin
from pwchem.constants import MGL_DIC
from autodock.objects import AutodockGrid

from schrodingerScipion.objects import SchrodingerGrid
from schrodingerScipion import Plugin as schrodinger_plugin

class ProtSchrodingerGridADT(EMProtocol):
    """Calls ADT to prepare a grid with an input that is a set of binding sites, a coordinate or a Schrodinger grid"""
    _label = 'grid definition using ADT'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct",
                       label='Input structure:', allowsNull=False)
        form.addParam('inputSites', MultiPointerParam, pointerClass="BindingSite, SetOfBindingSites, SchrodingerGrid",
                       label='Binding sites or grids:', allowsNull=False)
        form.addParam('inputLibrary', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Library of compounds:', allowsNull=False,
                       help='They must be either mol2 or pdbqt')
        form.addParam('stepSize', FloatParam, default=0.375, label='Step size (A)',
                      expertLevel=LEVEL_ADVANCED)

        line = form.addLine('Inner box (Angstroms)')
        line.addParam('innerX', IntParam, default=10, label='X')
        line.addParam('innerY', IntParam, default=10, label='Y')
        line.addParam('innerZ', IntParam, default=10, label='Z')

        line = form.addLine('Outer box (Angstroms)')
        line.addParam('outerX', IntParam, default=30, label='X')
        line.addParam('outerY', IntParam, default=30, label='Y')
        line.addParam('outerZ', IntParam, default=30, label='Z')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('preparationStep')
        # self._insertFunctionStep('createOutput')

    def getSiteCentroidSchrodingeBindingSite(self, fnSite):
        args = schrodinger_plugin.getPluginHome('utils/schrodingerUtils.py') + " centroid %s %s" % \
               (fnSite, self._getTmpPath("centroid.txt"))
        schrodinger_plugin.runSchrodinger(self, "python3", args)
        fh = open(self._getTmpPath("centroid.txt"))
        line = fh.readline()
        fh.close()
        x, y, z = line.split()
        return x,y,z

    def getGridCentroid(self, fnDir, fnLastDir):
        fnIn = os.path.join(fnDir, "%s.in" % fnLastDir)
        if not os.path.exists(fnIn):
            fnIn = os.path.join(fnDir, "%s.inp" % fnLastDir)
        fh = open(fnIn)
        for line in fh.readlines():
            if "GRID_CENTER" in line:
                tokens = line.replace(',', ' ').split()
                x = tokens[1].strip()
                y = tokens[2].strip()
                z = tokens[3].strip()
        fh.close()
        return x,y,z

    def prepareGrid(self, site, fnLigandBase):
        fnSite = site.getFileName()
        fnTarget = self.inputStructure.get().getFileName()
        if type(site)==BindingSite:
            n = fnSite.split('@')[0]
            x,y,z = self.getSiteCentroidSchrodingeBindingSite(fnSite)
        else:
            fnDir = os.path.split(fnSite)[0]
            fnLastDir = os.path.split(fnDir)[1]
            x, y, z = self.getGridCentroid(fnDir, fnLastDir)
            if "glide" in fnLastDir:
                n = fnLastDir.split("_")[1]
            else:
                n = fnLastDir.split("-")[1]

        fnGridDir = "grid-%s" % n
        fnGridDirAbs = self._getPath(fnGridDir)
        makePath(fnGridDirAbs)
        fnTargetLocal = os.path.join(fnGridDirAbs, os.path.split(fnTarget)[1])
        createLink(fnTarget, fnTargetLocal)

        fhTemplate = open(os.path.join(fnGridDirAbs,"inputParams.gpf"),'w')
        fhTemplate.write("npts %d %d %d\n"%(round(self.outerX.get()/self.stepSize.get()),
                                            round(self.outerY.get()/self.stepSize.get()),
                                            round(self.outerZ.get()/self.stepSize.get())))
        fhTemplate.write("spacing %f\n"%self.stepSize.get())
        fhTemplate.write("gridcenter %s %s %s\n"%(x,y,z))
        fhTemplate.close()

        args = " -r %s -l %s -d %s"%(fnTargetLocal,self._getExtraPath("Ligands/"+fnLigandBase),
                                     self._getExtraPath("Ligands"))
        libraryGPF = os.path.join(fnGridDirAbs,"library.gpf")
        args += " -o %s"%libraryGPF

        self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                    pwchem_plugin.getADTPath('Utilities24/prepare_gpf4.py')+args)

        args = "-p library.gpf -l library.glg"
        self.runJob(pwchem_plugin.getAutodockPath("autogrid4"),args, cwd=fnGridDirAbs)
        return fnGridDirAbs

    def preparationStep(self):
        # Create a link to all ligands
        fnLigandDir=self._getExtraPath("Ligands")
        makePath(fnLigandDir)
        for small in self.inputLibrary.get():
            fnSmall = small.getFileName()
            fnBase = os.path.split(fnSmall)[1]
            createLink(fnSmall,os.path.join(fnLigandDir,fnBase))

        # Now create the grids
        for site in self.inputSites:
            if type(site.get())==BindingSite or type(site.get())==SchrodingerGrid:
                gridDir = self.prepareGrid(site.get(), fnBase)
                self.createOutputSingle(site.get(), gridDir, site)
            elif type(site.get())==SetOfBindingSites:
                setOfSites = site
                for sitee in setOfSites.get():
                    gridDir = self.prepareGrid(sitee, fnBase)
                    self.createOutputSingle(sitee, gridDir, setOfSites)

    def createOutputSingle(self, site, fnDir, srcObj):
        score = None
        dscore = None
        if hasattr(site.get(), "score"):
            score = site.get().score.get()
            dscore = site.get().dscore.get()

        if os.path.exists(fnDir):
            gridFile = AutodockGrid(filename=fnDir)
            gridFile.structureFile = String(self.inputStructure.get().getFileName())
            if score:
                gridFile.bindingSiteScore = Float(score)
                gridFile.bindingSiteDScore = Float(dscore)

            n = fnDir.split('grid-')[1]
            outputDict = {'outputGrid%s' % n: gridFile}
            self._defineOutputs(**outputDict)
            self._defineSourceRelation(srcObj, gridFile)
