# **************************************************************************
# *
# * Authors: Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import STEPS_PARALLEL, PointerParam, BooleanParam, FloatParam, IntParam, LEVEL_ADVANCED
from pyworkflow.utils import redStr

# Plugin imports
from schrodingerScipion import Plugin

class ProtSchrodingerQikprop(EMProtocol):
	""" TODO: find out & place here """
	_label = 'qikprop'
	_program = ""

	# --------------------------- Class constructor --------------------------------------------
	def __init__(self, **args):
		# Calling parent class constructor
		super().__init__(**args)

		# Defining execution mode. Steps will take place in parallel now
		# Full tutorial on how to parallelize protocols can be read here:
		# https://scipion-em.github.io/docs/release-3.0.0/docs/developer/parallelization.html
		self.stepsExecutionMode = STEPS_PARALLEL

	def _defineParams(self, form):
		""" This function defines the params to be included in the protocol's form. """
		# Add parallel section
		form.addParallelSection(threads=4)

		# Create form
		form.addSection(label='Input')
		# Main params
		form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
			label='Input small molecules:', help='Input small molecules to be treated.')
		form.addParam('fast', BooleanParam, default=False, label='Run in fast mode:',
			help='Run job in fast mode (no dipole, homo or lumo).')
		form.addParam('sim', BooleanParam, default=True, label='Generate similar known drugs:',
			help='Generate a list of known drugs most similar to each processed molecule.')
		form.addParam('nsim', IntParam, default=5, label='Number of most similar drugs:', condition="sim==True",
			help="Number of most similar drug molecules to report.")
		
		# Additional optional params
		form.addParam('neut', BooleanParam, default=True, label='Neutralize molecules:', expertLevel=LEVEL_ADVANCED,
			help='Neutralize molecules in Maestro formatted files prior to processing.')
		form.addParam('altclass', BooleanParam, default=False, label='Alternative class:', expertLevel=LEVEL_ADVANCED,
			help='Run additional SASA, PSA calculations using an alternative class definition.')
		form.addParam('useAltprobe', BooleanParam, default=False, label='Use alternative probe:', expertLevel=LEVEL_ADVANCED,
			help='Run additional SASA, PSA calculations using an alternative probe radii (1.4A by default).')
		form.addParam('altprobe', FloatParam, default=1.4, label='Alternative probe:', expertLevel=LEVEL_ADVANCED, condition='useAltprobe==True')
		form.addParam('recap', BooleanParam, default=False, label='Replace CombiGlide with methyl:', expertLevel=LEVEL_ADVANCED,
			help='Replace the CombiGlide functional group with a methyl group prior to processing.\n'
				'When used in the CombiGlide reagent-preparation process, gives properties for the \'naked sidechain\'.')
		
		# Utils params
		form.addParam('cleanTmps', BooleanParam, default='True', label='Clean temporary files: ', expertLevel=LEVEL_ADVANCED,
			help='Clean temporary files after finishing the execution.\nThis is useful to reduce unnecessary disk usage.')

	# --------------------------- INSERT steps functions --------------------
	def _insertAllSteps(self):
		""" This function inserts all steps functions that will be run when running the protocol. """
		# Getting common command string for every call
		baseCommand = self.getQikpropBaseCmd()

		# For every SmallMolecule in the input set, build complete command
		for molecule in self.getInputFiles():
			self._insertFunctionStep('runQikpropStep', baseCommand, molecule)

	def runQikpropStep(self, baseCommand, molecule):
		""" This function runs the schrodinger binary file with the given params. """
		self.runJob(baseCommand, f' {molecule}', cwd=self._getExtraPath())

		# Clean tmp files if selected
		if self.cleanTmps.get():
			self.runJob('rm -rf', f' {self.getTmpFiles(molecule)}', cwd=self._getExtraPath())
	
	# --------------------------- INFO functions --------------------------------------------
	def _validate(self):
		"""
		The function of this hook is to add some validation before the protocol
		is launched to be executed. It should return a list of errors. If the list is
		empty the protocol can be executed.
		"""
		errors = []

		# Checking if MPI is selected (only threads are allowed)
		if self.numberOfMpi > 1:
			errors.append('MPI cannot be selected, because Scipion is going to drop support for it. Select threads instead.')
		
		# Checking if altprobe is used, and, if so, that the value is 0 or greater
		if self.useAltprobe.get() and self.altprobe.get() < 0:
			errors.append('Altprobe has to be a number equal or greater than 0.')

		return errors
	
	# --------------------------- Utils functions --------------------
	def getQikpropBaseCmd(self):
		""" This function returns the command string to run qikprop. """
		# Command starts with the executable file
		command = self.getQikpropBinaryFile()

		# Add permanent flags
		command += f' {self.getFastFlag()} {self.getSimFlag()} {self.getNeutFlag()} -LOCAL'

		# Add optional flags
		command += self.getNSim() + self.getAltClassFlag() + self.getAltProbeFlag() + self.getRecapFlag()

		# Return formatted command string
		return command

	def getQikpropBinaryFile(self):
		""" This function returns the location for the Schrodinger qikprop binary file. """
		# Getting path to the binary
		binaryPath = os.path.join(Plugin.getVar('SCHRODINGER_HOME'), 'qikprop')

		# If path exists, return it
		if os.path.exists(binaryPath):
			return binaryPath
		
		# If path was not found, raise exception
		raise FileNotFoundError(redStr(f"Path \"{binaryPath}\" not found. Is variable SCHRODINGER_HOME properly set within scipion.conf file?"))
	
	def getInputFiles(self):
		""" This function returns a list with the full path to each one of the input files. """
		return [os.path.abspath(molecule.getFileName()) for molecule in self.inputSmallMolecules.get()]
	
	def getFastFlag(self):
		""" This function returns the flag string corresponding to the fast flag param. """
		return '-fast' if self.fast.get() else '-nofast'

	def getSimFlag(self):
		""" This function returns the flag string corresponding to the sim flag param. """
		return '-sim' if self.sim.get() else '-nosim'
	
	def getNSim(self):
		""" This function returns the flag string corresponding to the nsim flag param. """
		return f' -nsim {self.nsim.get()}' if self.sim.get() else ''
	
	def getNeutFlag(self):
		""" This function returns the flag string corresponding to the neut flag param. """
		return '-neut' if self.neut.get() else '-noneut'
	
	def getAltClassFlag(self):
		""" This function returns the flag string corresponding to the altclass flag param. """
		return ' -altclass' if self.altclass.get() else ''
	
	def getAltProbeFlag(self):
		""" This function returns the flag string corresponding to the altprobe flag param. """
		return f' -altprobe {self.altprobe.get()}' if self.useAltprobe.get() else ''
	
	def getRecapFlag(self):
		""" This function returns the flag string corresponding to the sim recap param. """
		return ' -recap' if self.recap.get() else ''
	
	def getTmpFiles(self, molecule):
		""" This  function returns the temporary files related to the execution of qikprop for the given molecule. """
		print(molecule)
		return ''
