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
import os, shutil

# Scipion em imports
from pyworkflow.tests import setupTestProject, BaseTest
from pyworkflow.utils import yellowStr, redStr

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

		# Launching ligand preparation and obtaining output
		ligPrepProtocol = ligandPrep._runLigandPreparation()

		# Checking if raw small molecules were properly obtained
		rawSmallMols = getattr(ligPrepProtocol, 'inputSmallMolecules', None)
		if rawSmallMols is None:
			raise AssertionError(redStr("There was an error obtaining the input raw small molecules."))

		# Checking if ligand prepared small molecules were properly obtained
		ligprepSmallMols = getattr(ligPrepProtocol, LIGPREP_OUTPUTATTRIBUTE, None)
		if ligprepSmallMols is None:
			raise AssertionError(redStr("There was an error obtaining the input ligand prepared small molecules."))

		# Getting output details to input for import small molecules protocol
		# Checking if raw small molecules's path was properly obtained
		rawMoleculesPath = cls._getRawInputMolsPath(rawSmallMols)
		if rawMoleculesPath is None:
			raise AssertionError(redStr("There was an error obtaining the path for the input raw small molecules."))

		# Getting path to ligand prepared small molecules before path variables are overwritten by setupTestProject
		ligprepSmallMolsPath = os.path.abspath(ligPrepProtocol._getExtraPath())

		# Setting up project again to overwrite temp project variables
		setupTestProject(cls)

		# Importing small molecules
		cls.protRawImportSmallMols = cls._runImportSmallMols(rawMoleculesPath)
		cls.protLigPrepImportSmallMols = cls._runImportSmallMols(ligprepSmallMolsPath, processed=True)

	@classmethod
	def _runImportSmallMols(cls, moleculesPath, processed=False):
		# Creating and launching protocol
		protImportSmallMols = cls.newProtocol(
			ProtChemImportSmallMolecules,
			filesPath=moleculesPath,
			filesPattern=f'*.{"sdf" if processed else "mol2"}')
		cls.launchProtocol(protImportSmallMols, wait=True)
		return protImportSmallMols
	
	@classmethod
	def _runQikprop(cls, processed=False):
		""" This function creates, runs, and returns a Qikprop protocol. """
		inputProtocol = cls.protLigPrepImportSmallMols if processed else cls.protRawImportSmallMols
		protQikprop = cls.newProtocol(
			ProtSchrodingerQikprop,
			inputSmallMolecules=inputProtocol.outputSmallMolecules,
			fast=True,
			numberOfMpi=1,
			numberOfThreads=5
		)
		cls.launchProtocol(protQikprop, wait=True)
		return protQikprop
	
	@classmethod
	def _getRawInputMolsPath(cls, rawInputMols):
		""" This method returns the path for the given raw input small molecules. """
		try:
			return os.path.dirname(os.path.abspath(rawInputMols.get().getFirstItem().getFileName()))
		except AttributeError:
			return None
	
	@classmethod
	def _removeTmpElements(cls, tmpElements):
		""" This function removes all given temporary files and directories. """
		# Removing selected elements
		for item in tmpElements:
			if os.path.exists(item):
				if os.path.isdir(item):
					shutil.rmtree(item)
				else:
					os.remove(item)

	def test1(self):
		""" This function tests a qikprop execution when the proper input is received. """
		# Running Qikprop
		print(yellowStr("Running Qikprop with a ligand prepared set of small molecules."))
		qikpropProt = self._runQikprop(processed=True)
		self.assertIsNotNone(getattr(qikpropProt, QIKPROP_OUTPUTATTRIBUTE, None), msg="There was an error running Qikprop and it did not produce output.")
		
		# Getting summary (should be an empty list, as no warnings are produced)
		summaryString = '\n'.join(qikpropProt._summary())
		self.assertEqual(summaryString, '', msg=f"There was an error running Qikprop and some warnings were printed into summary:\n{summaryString}")

	def test2(self):
		""" This function tests a qikprop execution when no proper input is received. """
		# Running Qikprop
		print(yellowStr("Running Qikprop with a raw set of small molecules."))
		qikpropProt = self._runQikprop()
		self.assertIsNotNone(getattr(qikpropProt, QIKPROP_OUTPUTATTRIBUTE, None), msg="There was an error running Qikprop and it did not produce output.")
		
		# Getting summary (should be a list with an element for each warning produced)
		summaryString = '\n'.join(qikpropProt._summary())
		self.assertNotEqual(summaryString, '', msg="There was an error running Qikprop and no warnings were printed into summary when there should be some.")

		# Last test calls cleaning function so it does not count as a separate test
		self._removeTmpElements(self.tmpElements)