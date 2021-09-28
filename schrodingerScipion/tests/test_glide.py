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

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1_noHETATM.pdb'))
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runImportSmallMols(cls):
      protImportSmallMols = cls.newProtocol(
        ProtChemImportSmallMolecules,
        filesPath=cls.dsLig.getFile('mol2'))
      cls.launchProtocol(protImportSmallMols)
      cls.protImportSmallMols = protImportSmallMols

    def _runTargetPreparation(self, kwargs):
        protPrepWizard = self.newProtocol(
            ProtSchrodingerPrepWizard,
            inputStructure=self.protImportPDB.outputPdb,
            **kwargs)
        self.launchProtocol(protPrepWizard)
        return protPrepWizard

    def _runLigandPreparation(self):
        protPrepLigand = self.newProtocol(
            ProtSchrodingerLigPrep,
            inputSmallMols=self.protImportSmallMols.outputSmallMols,
            ionization=1)
        self.launchProtocol(protPrepLigand)
        return protPrepLigand

    def _runSitemap(self, targetProt):
        protSitemap = self.newProtocol(
            ProtSchrodingerSiteMap,
            inputStructure=targetProt.outputStructure)

        self.launchProtocol(protSitemap)
        pocketsOut = getattr(protSitemap, 'outputPockets', None)
        self.assertIsNotNone(pocketsOut)
        return protSitemap
    
    def _runFilterSites(self, siteProt):
        protFilter = self.newProtocol(
            ProtSetFilter,
            operation=ProtSetFilter.CHOICE_RANKED,
            threshold=2,
            rankingField='_score')
        protFilter.inputSet.set(siteProt)
        protFilter.inputSet.setExtended('outputPockets')

        self.launchProtocol(protFilter)
        return protFilter

    def _runGridDefinition(self, filterProt, targetProt):
        protGrid = self.newProtocol(
            ProtSchrodingerGridSiteMap,
            inputSchAtomStruct=targetProt.outputStructure,
            innerAction=1, diameterNin=1.2,
            outerAction=1, diameterNout=0.8)
        protGrid.inputSetOfPockets.set(filterProt)
        protGrid.inputSetOfPockets.setExtended("outputPockets")

        self.launchProtocol(protGrid)
        gridsOut = getattr(protGrid, 'outputGrids', None)
        self.assertIsNotNone(gridsOut)
        return protGrid

    def _runGlideDocking(self, ligProt, gridProt):
        protGlide = self.newProtocol(
            ProtSchrodingerGlideDocking,
            inputGridSet=gridProt.outputGrids,
            inputLibrary=ligProt.outputSmallMols,
            doConvertOutput=True, convertType=1)

        self.launchProtocol(protGlide)
        return protGlide

    def testGlide(self):
        prepProt = self._runTargetPreparation(self.getPrepTargetWizardArgs())
        siteProt = self._runSitemap(prepProt)
        filterProt = self._runFilterSites(siteProt)
        gridProt = self._runGridDefinition(filterProt, prepProt)

        ligProt = self._runLigandPreparation()
        glideProt = self._runGlideDocking(ligProt, gridProt)

    def getPrepTargetWizardArgs(self):
        args={'stage1':True, 'hydrogens': 1,
              'stage2':True, 'propKa':True,
              'stage3':True,
              'stage4':True, 'ms':True}
        return args




