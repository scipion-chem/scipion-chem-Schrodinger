# **************************************************************************
# *
# * Authors:		 Martín Salinas (martin.salinas@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
# General imports
import os

# Scipion em imports
from pyworkflow.tests import setupTestProject, BaseTest
from pyworkflow.utils import yellowStr

# Scipion chem imports
from pwchem.protocols import ProtChemImportSmallMolecules

# Plugin imports
from ..protocols import ProtSchrodingerQikprop
from ..protocols.protocol_ligprep import OUTPUTATTRIBUTE as LIGPREP_OUTPUTATTRIBUTE
from ..protocols.protocol_qikprop import OUTPUTATTRIBUTE as QIKPROP_OUTPUTATTRIBUTE
from ..tests import TestSchroLigPrep

class TestSchroQikprop(BaseTest):
	""" This test checks if qikprop protocol works correctly. """
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)
		cls.tmpElements = []

		# Getting ProtSchrodingerLigPrep test so it can produce necessary inputs for this test
		ligandPrep = TestSchroLigPrep()
		ligandPrepPath = cls.getOutputPath().replace(cls.__name__, type(ligandPrep).__name__)

		# Check if temp project already exists
		if not os.path.isdir(ligandPrepPath):
			# If not, add to delete after end of this test so no extra files and directories are left behind
			cls.tmpElements.append(ligandPrepPath)
		
		# Launching TestSchroLigPrep and obtaining output
		# running import small molecules and obtainig input set
		print(yellowStr("Launching Ligand preparation test to reuse its outputs as inputs for qikprop."))
		ligandPrep.setUpClass()
		#rawSmallMols = getattr(ligandPrep.protImportSmallMols, 'outputSmallMolecules', None)
		#cls.assertIsNotNone(rawSmallMols, "There was an error obtaining the raw small molecules.")
		# Launching ligand preparation and obtaining output
		ligPrepProtocol = ligandPrep._runLigandPreparation()
		ligprepSmallMols = getattr(ligPrepProtocol, LIGPREP_OUTPUTATTRIBUTE, None)
		cls.assertIsNotNone(ligprepSmallMols, "There was an error obtaining the ligand prepared small molecules.")

		# Getting output details to input for import small molecules protocol
		#rawMoleculesPath = os.path.abspath(ligPrepProtocol.inputSmallMolecules.get())
		ligPrepMoleculesPath = os.path.abspath(ligPrepProtocol._getExtraPath())

		# Setting up project again to overwrite temp project variables
		setupTestProject(cls)

		# Importing small molecules
		#cls.rawMols = cls._runImportSmallMols(cls, rawMoleculesPath)
		cls.ligPrepMols = cls._runImportSmallMols(ligPrepMoleculesPath, processed=True)

	@classmethod
	def _runImportSmallMols(cls, moleculesPath, processed=False):
		cls.protImportSmallMols = cls.newProtocol(
			ProtChemImportSmallMolecules,
			filesPath=moleculesPath,
			filesPattern=f'*.{"sdf" if processed else "mol2"}')
		cls.launchProtocol(cls.protImportSmallMols)
	
	@classmethod
	def _runQikprop(cls):
		""" This function creates, runs, and returns a Qikprop protocol. """
		protQikprop = cls.newProtocol(
			ProtSchrodingerQikprop,
			inputSmallMolecules=cls.protImportSmallMols.outputSmallMolecules,
			fast=True)
		cls.launchProtocol(protQikprop)
		return protQikprop

	def test1(self):
		""" This function tests a qikprop execution with the proper input received. """
		print(yellowStr("Running Qikprop with a ligand prepared set of small molecules."))
		qikpropProt = self._runQikprop()
		self.assertIsNotNone(getattr(qikpropProt, QIKPROP_OUTPUTATTRIBUTE, None))
