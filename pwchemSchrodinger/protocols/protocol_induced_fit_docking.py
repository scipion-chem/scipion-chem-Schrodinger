# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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
# General imports
import random as rd
import glob, os

# Scipion em imports
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import Group, EnumParam, PointerParam, StringParam, FloatParam, LabelParam, TextParam, IntParam, LEVEL_ADVANCED
from pyworkflow.utils import Message

# Scipion chem imports
from pwchem.utils import natural_sort

# Plugin imports
from .. import Plugin as schrodingerPlugin
from ..constants import MSJ_SYSMD_INIT
from ..protocols.protocol_glide_docking import ProtSchrodingerGlideDocking
from ..utils import createMSJDic, buildSimulateStr, setAborted

jobControlProg = schrodingerPlugin.getHome('jobcontrol')
structConvertProg = schrodingerPlugin.getHome('utilities/structconvert')
maeSubsetProg = schrodingerPlugin.getHome('utilities/maesubset')
mmgbsaProg = schrodingerPlugin.getHome('prime_mmgbsa')

#  0       1     2       3      4      5      6       7     8     9      10     11    12    13    14
COMPILE, GLIDE, IDOCK, PPREP, PFLEX, PENER, PHELIX, PLOOP, PMIN, PREF, PSIDE, SCORE, SORT, TRIM, VDW = list(range(15))

class ProtSchrodingerIFD(ProtSchrodingerGlideDocking):
  """Perform an Induced Fit Docking protocol using the user defined set of stages"""
  _label = 'Induced Fit Docking (IFD)'

  stageTypes = {'Compile residues': 'COMPILE_RESIDUE_LIST',
                'Glide': 'GLIDE_DOCKING2', 'Initial docking': 'INITIAL_DOCKING', 'Receptor preparation': 'PPREP',
                'Predict flexibility': 'PREDICT_FLEXIBILITY', 'Energy calculation': 'PRIME_ECALC',
                'Helix refinement': 'PRIME_HELIX', 'Loop refinement': 'PRIME_LOOP',
                'Residues minimization': 'PRIME_MINIMIZATION', 'Residues refinement': 'PRIME_REFINEMENT',
                'Side-chain prediction': 'PRIME_SIDECHAIN',
                'Score poses': 'SCORING', 'Sort and filter poses': 'SORT_AND_FILTER',
                'Trim side-chains': 'TRIM_SIDECHAINS', 'Estimate side-chain flexibility': 'VDW_SCALING'}

  _omitParamNames = ['runName', 'runMode', 'insertStep', 'summarySteps', 'deleteStep', 'watchStep',
                     'workFlowSteps', 'hostName', 'numberOfThreads', 'numberOfMpi',
                     'inputSmallMolecules', 'fromPockets', 'inputAtomStruct', 'radius', 'inputStructROIs',
                     'inputGridSet']

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def getStageParamsDic(self, type='All'):
    '''Return a dictionary as {paramName: param} of the stage parameters of the formulary.
    Type'''
    paramsDic = {}
    for paramName, param in self._definition.iterAllParams():
      if paramName not in self._omitParamNames and not isinstance(param, Group):
        if type == 'All' or (type == 'Enum' and isinstance(param, EnumParam)) or (type == 'Normal' and not isinstance(param, EnumParam)):
          paramsDic[paramName] = param
    return paramsDic

  def _defineParams(self, form):
    # Defining costant variables
    residueListHelp = 'Comma-separated list of residues to add to the refinement list. These should be residues that lie outside the distance cutoff'
    wizardHelp = 'Click on the wizard to save the specified residues in the list below'
    cutoffDistanceLabel = 'Cutoff distance:'

    form.addSection(label=Message.LABEL_INPUT)
    form = self._defineGlideReceptorParams(form)
    group = form.addGroup('Ligands')
    group.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                   label='Input small molecules:', help='Input small molecules to be docked with IFD')

    form.addSection('IFD stages')
    group = form.addGroup('Add Stage')
    group.addParam('stageType', EnumParam, default=0, label='Stage type: ', choices=list(self.stageTypes.keys()),
                   help='Type of the next stage to add to the protocol. For details about the different stages, '
                        'check the Schrodinger docs in {}'.
                   format(schrodingerPlugin.getHome('docs/inducedfit_command_reference/ifd_command_infile.htm')))

    # COMPILE_RESIDUE_LIST stage parameters
    group.addParam('residuesDefinition', EnumParam, label='Define residues by: ', default=1,
                   choices=['Distance form center', 'Residue list'], condition=f'stageType == {COMPILE}',
                   help='Whether to determine the list of residues to compile by using a center and a cutoff distance '
                        'or directly a list of residues')
    group.addParam('residueCenter', StringParam, default='Z:999', label='Center to determine residues: ',
                   condition=f'stageType == {COMPILE} and residuesDefinition == 0',
                   help='List of residues from which to measure the cutoff distance. Default: Z:999, which is the '
                        'default for the ligand.')
    group.addParam('residuesCutoff', FloatParam, label=cutoffDistanceLabel,
                   default=5.0, condition=f'stageType == {COMPILE} and residuesDefinition == 0',
                   help='Cutoff distance (in angstroms) from the ligand pose, within which residues that have any '
                        'atoms are included in the refinement list')

    group.addParam('addResAdd', LabelParam, label='Insert residues to add list: ',
                   condition=f'stageType == {COMPILE} and residuesDefinition == 1',
                   help=wizardHelp)
    group.addParam('residuesAdd', TextParam, default='', width=120,
                   label='Compile residues to add: ', condition=f'stageType == {COMPILE} and residuesDefinition == 1',
                   help=residueListHelp)
    group.addParam('addResOmit', LabelParam, label='Insert residues to omit list: ',
                   condition=f'stageType == {COMPILE} and residuesDefinition == 1',
                   help=wizardHelp)
    group.addParam('residuesOmit', TextParam, default='', width=120,
                   label='Compile residues to omit: ', condition=f'stageType == {COMPILE} and residuesDefinition == 1',
                   help='Comma-separated list of residues to omit to the refinement list')

    # GLIDE_DOCKING2 and INITIAL_DOCKING stage parameters
    form = self._defineGlideParams(form, condition=f'stageType in [{GLIDE}, {IDOCK}]')

    # PPREP stage parameters
    group.addParam('convergenceRMSD', FloatParam, label='Convergence RMSD: ',
                   default=0.2, condition=f'stageType == {PPREP}',
                   help='RMSD value which specifies the convergence threshold for the constrained minimization of '
                        'the Glide protein preparation. The minimization ends when the RMSD is less than or equal to '
                        'value.')

    # PRIME_HELIX and PRIME_LOOP stage parameters
    group.addParam('residuesRegion', StringParam, default='',
                   label='Residues to take into account: ', condition=f'stageType in [{PHELIX}, {PLOOP}]',
                   help='Residues of mobile region that contains the helix or loop')
    group.addParam('residuesHelix', StringParam, default='',
                   label='Helix residues: ', condition=f'stageType in [{PHELIX}]',
                   help='Residue of the helix')
    group.addParam('helixLoopCutoff', FloatParam, label=cutoffDistanceLabel,
                   default=5.0, condition=f'stageType in [{PHELIX}, {PLOOP}]',
                   help='HELIX: Distance cutoff for prediction of side chains close to helix.\n'
                        'LOOP: Threshold for inclusion of residues for side-chain refinement. Any residues with atoms '
                        'within this distance are included')
    group.addParam('maxEnergyGap', FloatParam, label=cutoffDistanceLabel,
                   default=10000.0, condition=f'stageType in [{PHELIX}, {PLOOP}]',
                   help='Maximum energy gap for saved structures.'
                        'Energy threshold for predicted loop structures (in kcal/mol). Structures are discarded if '
                        'their energy is more than this amount above the lowest-energy structure.')
    group.addParam('maxStructures', IntParam, label='Max. number of structures: ',
                   default=1000, condition=f'stageType in [{PHELIX}, {PLOOP}]',
                   help='Maximum number of structures to store')

    # PRIME_REFINEMENT and PRIME_MINIMIZATION stage parameters
    group.addParam('nMinPasses', IntParam, label='Number of refinement passes: ',
                   default=1, condition=f'stageType in [{PREF}]',
                   help='In Prime refinement, the side chains of the residue list compiled previously are optimized, '
                        'then the residues are minimized. This parameter controls the number of times to repeat those '
                        'two steps')

    # SCORING stage parameters
    group.addParam('scoreName', StringParam, default='',
                   label='Name of property to add to Maestro files: ', condition=f'stageType in [{SCORE}]',
                   help='Name of property to add to Maestro files')
    group.addParam('scoreTerms', TextParam, label='Scoring terms: ', default='', width=120,
                   condition=f'stageType in [{SCORE}]',
                   help='Add terms to the scoring function. You can include multiple TERM keywords to define the '
                        'scoring function (one per line). Format: coeff,property,stage, where coeff is the '
                        'coefficient, property is the property from the Maestro output file, and stage is the index '
                        'of the property-generating stage, counting backwards in the input file with 0 for the '
                        'previous stage. \nSee for details: '
                        '{schrodingerPlugin.getHome("docs/inducedfit_command_reference/ifd_command_infile_scoring.htm")}')

    # SORT_AND_FILTER stage parameters
    group.addParam('poseFilter', StringParam, default='',
                   label='Property to filter poses: ', condition=f'stageType in [{SORT}]',
                   help='Name of Maestro property for filtering poses, for example, r_psp_Prime_Energy')
    group.addParam('poseKeep', StringParam, default='',
                   label='Threshold on property for filtering poses: ', condition=f'stageType in [{SORT}]',
                   help='Threshold on property for filtering poses. The syntax is as follows: '
                        'n% Keep the n% of poses with the lowest property values'
                        'n# Keep the n poses with the lowest property values'
                        'n Keep poses with property values within n of the lowest value.')
    group.addParam('ligandFilter', StringParam, default='',
                   label='Property to filter ligands: ', condition=f'stageType in [{SORT}]',
                   help='Name of Maestro property for filtering ligands, for example, r_psp_Prime_Energy')
    group.addParam('ligandKeep', StringParam, default='',
                   label='Threshold on property for filtering ligands: ', condition=f'stageType in [{SORT}]',
                   help='Threshold on property for filtering ligands. The syntax is as follows: '
                        'n% Keep the n% of poses with the lowest property values'
                        'n# Keep the n poses with the lowest property values'
                        'n Keep poses with property values within n of the lowest value.')

    # TRIM_SIDECHAINS stage parameters
    group.addParam('trimResidues', EnumParam, label='Residues to mutate to Ala: ', default=1,
                   display=EnumParam.DISPLAY_HLIST,
                   choices=['Manual specification', 'Automatic'], condition=f'stageType == {TRIM}',
                   help='Residues for temporary mutation to alanine. Either a list of residues must be given, or AUTO.')
    group.addParam('trimMethod', EnumParam, label='Method for automatic specification of residues: ', default=0,
                   choices=['BFACTOR', 'FLEXIBILITY'], condition=f'stageType == {TRIM} and trimResidues == 1',
                   help='Method for automatic specification of residues. '
                        'BFACTOR: Mutate residues near the ligand with the highest B factors, subject to a cutoff on '
                        'the B-factor and the number of residues.\n'
                        'FLEXIBILITY: Mutate residues for which other rotamers exist that don’t have steric '
                        'clashes with the rest of the protein.')
    group.addParam('cutOffBFactor', FloatParam, label='B-factor cutoff: ',
                   default=40.0, condition=f'stageType == {TRIM} and trimResidues == 1 and trimMethod == 0',
                   help='B-factor cutoff above which residues are selected for mutation.')
    group.addParam('maxBFactor', IntParam, label='Maximum number of residues to mutate: ',
                   default=3, condition=f'stageType == {TRIM} and trimResidues == 1 and trimMethod == 0',
                   help='Maximum number of residues to mutate. If more than N residues have B-factors greater than the '
                        'value of BFACTOR_CUTOFF, the N residues with the highest B-factors are selected')
    group.addParam('maxFlexibility', FloatParam, label='Max flexibility cutoff (Å): ',
                   default=5.0, condition=f'stageType == {TRIM} and trimResidues == 1 and trimMethod == 1',
                   help='If any heavy atom in the side chain is capable of moving more than this distance, the side '
                        'chain is considered flexible and is mutated')
    group.addParam('resolFlexibility', FloatParam, label='Granularity of rotamer search: ',
                   default=30.0, condition=f'stageType == {TRIM} and trimResidues == 1 and trimMethod == 1',
                   help='Granularity of the rotamer search for steric clashes, in degrees.')

    group.addParam('addResTrim', LabelParam, label='Insert residues to trim list: ',
                   condition=f'stageType == {TRIM} and trimResidues == 0',
                   help=wizardHelp)
    group.addParam('residuesTrim', TextParam, default='', width=120,
                   label='Residues to trim: ', condition=f'stageType == {TRIM} and trimResidues == 0',
                   help='List of residues to trim')

    group.addParam('extraParams', TextParam, label='Extra parameters: ', default='', width=120,
                   expertLevel=LEVEL_ADVANCED,
                   help='Any extra parameters to add in the current stage. Add each parameter in a different line as:\n'
                        'PARAMETER_NAME PARAMETER_VALUE')

    group = form.addGroup('Residue selection', condition=f'(stageType == {COMPILE} and residuesDefinition == 1) or (stageType in [{PHELIX}, {PLOOP}]) or (stageType == {TRIM} and trimResidues == 0)')
    group.addParam('selChain', StringParam, default='',
                   label='Chain of residues: ',
                   condition=f'(stageType == {COMPILE} and residuesDefinition == 1) or (stageType in [{PHELIX}, {PLOOP}]) or (stageType == {TRIM} and trimResidues == 0)',
                   help=residueListHelp)
    group.addParam('selResidue', StringParam, default='',
                   label='Chain of residues: ',
                   condition=f'(stageType == {COMPILE} and residuesDefinition == 1) or (stageType in [{PHELIX}, {PLOOP}]) or (stageType == {TRIM} and trimResidues == 0)',
                   help=residueListHelp)

    group = form.addGroup('Summary')
    group.addParam('insertStep', StringParam, default='',
                   label='Insert stage number: ',
                   help='Insert the defined protocol stage into the workflow on the defined position.\n'
                        'The default (when empty) is the last position')
    group.addParam('summarySteps', TextParam, width=120, readOnly=True,
                   label='Summary of stages', help='Summary of the defined stages. \nManual modification will have no '
                                                   'effect, use the wizards to add / delete the stages')
    group.addParam('deleteStep', StringParam, default='',
                   label='Delete stage number: ',
                   help='Delete the stage of the specified index from the workflow.')
    group.addParam('watchStep', StringParam, default='',
                   label='Watch stage number: ',
                   help='''Watch the parameters stage of the specified index from the workflow..\n
                       This might be useful if you want to change some parameters of a predefined stage.\n
                       However, the parameters are not changed until you add the new stage (and probably\n
                       you may want to delete the previous unchanged stage)''')
    group.addParam('workFlowSteps', StringParam, label='User transparent', condition='False')

    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
    # self._insertFunctionStep('simulationStep')
    # self._insertFunctionStep('createOutputStep')
    pass


  #######################################

  def addDefaultForMissing(self, msjDic):
    '''Add default values for missing parameters in the msjDic'''
    paramDic = self.getStageParamsDic()
    for pName in paramDic.keys():
      if pName not in msjDic:
        msjDic[pName] = paramDic[pName].default
    return msjDic

  ############# UTILS
  def buildMSJ_str(self):
    # todo: write the file neccessary for IFD
    '''Build the .msj (file used by IFD to specify the jobs performed by Schrodinger)
        defining the input parameters'''
    msjStr = MSJ_SYSMD_INIT

    if self.workFlowSteps.get() in ['', None]:
      msjDic = createMSJDic(self)
      msjStr += buildSimulateStr(self, msjDic)
    else:
      workSteps = self.workFlowSteps.get().split('\n')
      if '' in workSteps:
        workSteps.remove('')
      for wStep in workSteps:
        msjDic = eval(wStep)
        msjStr += buildSimulateStr(self, msjDic)

    return msjStr

  def createGUISummary(self):
    with open(self._getExtraPath("summary.txt"), 'w') as f:
      if self.workFlowSteps.get():
        f.write(self.createSummary())
      else:
        f.write(self.createSummary(createMSJDic(self)))

  def createSummary(self, msjDic=None):
    '''Creates the displayed summary from the internal state of the steps'''
    sumStr = ''
    if not msjDic:
      for i, dicLine in enumerate(self.workFlowSteps.get().split('\n')):
        if dicLine != '':
          msjDic = eval(dicLine)
          msjDic = self.addDefaultForMissing(msjDic)
          sumStr += f'{i+1}) {msjDic["stageType"]}'
          sumStr += '\n'
    else:
      msjDic = self.addDefaultForMissing(msjDic)
      sumStr += f'{msjDic["stageType"]}'
      sumStr += '\n'
    return sumStr

  def countSteps(self):
    stepsStr = self.summarySteps.get() if self.summarySteps.get() is not None else ''
    steps = stepsStr.split('\n')
    return len(steps) - 1

  def setAborted(self):
    super().setAborted()
    setAborted(self, jobControlProg)
