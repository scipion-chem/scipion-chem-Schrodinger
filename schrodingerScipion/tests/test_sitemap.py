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
from pwem.protocols import ProtImportPdb
from ..protocols import ProtSchrodingerSiteMap, ProtSchrodingerPrepWizard

class TestSitemap(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)
        cls._runImportPDB()

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1_noHETATM.pdb'))
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    def _runTargetPreparation(self, kwargs):
        protPrepWizard = self.newProtocol(
            ProtSchrodingerPrepWizard,
            inputAtomStruct=self.protImportPDB.outputPdb,
            **kwargs)
        self.launchProtocol(protPrepWizard)
        return protPrepWizard

    def _runSitemap(self, targetProt):
        protSitemap = self.newProtocol(
            ProtSchrodingerSiteMap,
            inputStructure=targetProt.outputStructure)

        self.launchProtocol(protSitemap)
        pdbOut = getattr(protSitemap, 'outputStructROIs', None)
        self.assertIsNotNone(pdbOut)

    def testSitemap(self):
        prepArgs = self.getPrepWizardArgs()
        prepProt = self._runTargetPreparation(prepArgs)
        self._runSitemap(prepProt)

    def getPrepWizardArgs(self):
        args={'stage1':True, 'hydrogens': 1,
              'stage2':True, 'propKa':True,
              'stage3':True,
              'stage4':True, 'ms':True}
        return args




