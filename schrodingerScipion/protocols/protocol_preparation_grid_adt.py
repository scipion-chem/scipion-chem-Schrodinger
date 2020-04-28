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

from pyworkflow.protocol.params import MultiPointerParam, BooleanParam, FloatParam, IntParam, PointerParam
from pyworkflow.object import String, Float
from pyworkflow.utils.path import createLink, makePath
from pyworkflow.protocol.constants import LEVEL_ADVANCED

from pwem.protocols import EMProtocol
from schrodingerScipion.objects import SchrodingerGrid
from bioinformatics.objects import BindingSite, SetOfBindingSites
from bioinformatics import Plugin as bioinformatics_plugin
from schrodingerScipion import Plugin as schrodinger_plugin

class ProtSchrodingerGridADT(EMProtocol):
    """Calls ADT to prepare a grid with an input that is a set of binding sites, a coordinate or a Schrodinger grid"""
    _label = 'grid definition using ADT'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('useCoordinates', BooleanParam, default=False,
                      label='Use coordinates')
        line = form.addLine('Center (Angstroms)', condition='useCoordinates')
        line.addParam('X', IntParam, default=0, label='X')
        line.addParam('Y', IntParam, default=0, label='Y')
        line.addParam('Z', IntParam, default=0, label='Z')

        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct",
                       label='Input structure:', allowsNull=False)
        form.addParam('inputSites', MultiPointerParam, pointerClass="BindingSite, SetOfBindingSites, SchrodingerGrid",
                       condition='not useCoordinates', label='Binding sites or grids:', allowsNull=False)
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
        if self.useCoordinates:
            x = str(self.X.get())
            y = str(self.Y.get())
            z = str(self.Z.get())
        if type(site)==BindingSite:
            n = fnSite.split('@')[0]
            if not self.useCoordinates:
                x,y,z = self.getSiteCentroidSchrodingeBindingSite(fnSite)
        else:
            fnDir = os.path.split(fnSite)[0]
            fnLastDir = os.path.split(fnDir)[1]
            if not self.useCoordinates:
                x, y, z = self.getGridCentroid(fnDir, fnLastDir)
            if "glide" in fnLastDir:
                n = fnLastDir.split("_")[1]
            else:
                n = fnLastDir.split("-")[1]
            print(fnDir)
            print(fnTarget)

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
        args += " -o %s"%os.path.join(fnGridDirAbs,"library.gpf")

        self.runJob(bioinformatics_plugin.getMGLPath('bin/pythonsh'),
                    bioinformatics_plugin.getADTPath('Utilities24/prepare_gpf4.py')+args)

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
                self.prepareGrid(site.get(), fnBase)
            elif type(site.get())==SetOfBindingSites:
                setOfSites = site
                for sitee in setOfSites.get():
                    self.prepareGrid(sitee, fnBase)

    def createOutputSingle(self, fnSite, fnStructureFile, score, dscore, srcObj):
        n = fnSite.split('@')[0]

        fnDir = self._getPath("grid-%s" % n)
        if os.path.exists(fnDir):
            fnBase = os.path.split(fnDir)[1]
            fnGrid = os.path.join(fnDir, '%s.zip' % fnBase)
            if os.path.exists(fnGrid):
                gridFile = SchrodingerGrid(filename=fnGrid)
                gridFile.structureFile = String(fnStructureFile)
                gridFile.bindingSiteScore = Float(score)
                gridFile.bindingSiteDScore = Float(dscore)

                n = fnDir.split('grid-')[1]
                outputDict = {'outputGrid%s' % n: gridFile}
                self._defineOutputs(**outputDict)
                self._defineSourceRelation(srcObj, gridFile)

    def createOutput(self):
        for site in self.inputSites:
            if type(site)==BindingSite:
                fnSite = site.get().getFileName()
                fnStructureFile = site.get().structureFile.get()
                score = site.get().score.get()
                dscore = site.get().dscore.get()
                self.createOutputSingle(fnSite, fnStructureFile, score, dscore, site)
            else:
                setOfSites = site
                for sitee in setOfSites.get():
                    fnSite = sitee.getFileName()
                    fnStructureFile = sitee.structureFile.get()
                    score = sitee.score.get()
                    dscore = sitee.dscore.get()
                    self.createOutputSingle(fnSite, fnStructureFile, score, dscore, setOfSites)
