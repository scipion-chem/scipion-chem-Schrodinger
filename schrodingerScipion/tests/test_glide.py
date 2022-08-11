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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb, ProtSetFilter
from pwchem.protocols import ProtChemImportSmallMolecules
from ..protocols import ProtSchrodingerSiteMap, ProtSchrodingerPrepWizard, \
    ProtSchrodingerLigPrep, ProtSchrodingerGridSiteMap, ProtSchrodingerGlideDocking

class TestGlideDocking(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet("smallMolecules")

        setupTestProject(cls)
        cls._runImportPDB()
        cls._runImportSmallMols()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)
        cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)

        cls.prepProt = cls._runTargetPreparation(getPrepTargetWizardArgs())
        cls.ligProt = cls._runLigandPreparation()
        cls._waitOutput(cls.prepProt, 'outputStructure', sleepTime=5)
        
    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0,
            pdbId='4erf')
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runImportSmallMols(cls):
      cls.protImportSmallMols = cls.newProtocol(
        ProtChemImportSmallMolecules,
        filesPath=cls.dsLig.getFile('mol2'))
      cls.launchProtocol(cls.protImportSmallMols)
      cls.proj.launchProtocol(cls.protImportSmallMols, wait=False)

    @classmethod
    def _runTargetPreparation(cls, kwargs):
        protPrepWizard = cls.newProtocol(
            ProtSchrodingerPrepWizard,
            cleanPDB=True, waters=False, rchains=True,
            chain_name='{"model": 0, "chain": "C", "residues": 93}',
            **kwargs)
        protPrepWizard.inputAtomStruct.set(cls.protImportPDB)
        protPrepWizard.inputAtomStruct.setExtended('outputPdb')
        cls.proj.launchProtocol(protPrepWizard, wait=False)
        return protPrepWizard

    @classmethod
    def _runLigandPreparation(cls):
        protPrepLigand = cls.newProtocol(
            ProtSchrodingerLigPrep,
            inputSmallMols=cls.protImportSmallMols.outputSmallMolecules,
            ionization=1)
        cls.proj.launchProtocol(protPrepLigand, wait=False)
        return protPrepLigand

    @classmethod
    def _runSitemap(cls, targetProt):
        protSitemap = cls.newProtocol(
            ProtSchrodingerSiteMap,
            inputStructure=targetProt.outputStructure)

        cls.launchProtocol(protSitemap)
        return protSitemap

    @classmethod
    def _runGridDefinition(cls, filterProt):
        protGrid = cls.newProtocol(
            ProtSchrodingerGridSiteMap,
            innerAction=1, diameterNin=0.8,
            outerAction=1, diameterNout=1.2)
        protGrid.inputStructROIs.set(filterProt)
        protGrid.inputStructROIs.setExtended("outputStructROIs")

        cls.launchProtocol(protGrid)
        return protGrid

    def _runGlideDocking(self, ligProt, gridProt):
        protGlide = self.newProtocol(
            ProtSchrodingerGlideDocking,
            doConvertOutput=True, convertType=1)

        protGlide.inputGridSet.set(gridProt)
        protGlide.inputGridSet.setExtended('outputGrids')
        protGlide.inputLibrary.set(ligProt)
        protGlide.inputLibrary.setExtended('outputSmallMolecules')

        self.launchProtocol(protGlide)
        return protGlide

    def testGlide(self):
        siteProt = self._runSitemap(self.prepProt)
        gridProt = self._runGridDefinition(siteProt)

        self._waitOutput(self.ligProt, 'outputSmallMolecules', sleepTime=5)
        glideProt = self._runGlideDocking(self.ligProt, gridProt)


def getPrepTargetWizardArgs():
    args={'stage1':True, 'hydrogens': 1,
          'stage2':True, 'propKa':True,
          'stage3':True,
          'stage4':True, 'ms':True}
    return args




