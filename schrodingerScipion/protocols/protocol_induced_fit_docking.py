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

import random as rd
import glob, os
from subprocess import check_call

from pyworkflow.protocol.params import *
from pwem.protocols import EMProtocol

from pwchem.utils import natural_sort

from .. import Plugin as schrodingerPlugin
from ..constants import *
from ..objects import SchrodingerSystem
from ..protocols.protocol_glide_docking import ProtSchrodingerGlideDocking

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
      if not paramName in self._omitParamNames and not isinstance(param, Group):
        if type == 'All':
          paramsDic[paramName] = param
        elif type == 'Enum' and isinstance(param, EnumParam):
          paramsDic[paramName] = param
        elif type == 'Normal' and not isinstance(param, EnumParam):
          paramsDic[paramName] = param
    return paramsDic

  def _defineParams(self, form):
    form.addSection(label='Input')
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
    group.addParam('residuesCutoff', FloatParam, label='Cutoff distance: ',
                   default=5.0, condition=f'stageType == {COMPILE} and residuesDefinition == 0',
                   help='Cutoff distance (in angstroms) from the ligand pose, within which residues that have any '
                        'atoms are included in the refinement list')

    group.addParam('addResAdd', LabelParam, label='Insert residues to add list: ',
                   condition=f'stageType == {COMPILE} and residuesDefinition == 1',
                   help='Click on the wizard to save the specified residues in the list below')
    group.addParam('residuesAdd', TextParam, default='', width=120,
                   label='Compile residues to add: ', condition=f'stageType == {COMPILE} and residuesDefinition == 1',
                   help='Comma-separated list of residues to add to the refinement list. These should be residues '
                        'that lie outside the distance cutoff')
    group.addParam('addResOmit', LabelParam, label='Insert residues to omit list: ',
                   condition=f'stageType == {COMPILE} and residuesDefinition == 1',
                   help='Click on the wizard to save the specified residues in the list below')
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
    group.addParam('helixLoopCutoff', FloatParam, label='Cutoff distance: ',
                   default=5.0, condition=f'stageType in [{PHELIX}, {PLOOP}]',
                   help='HELIX: Distance cutoff for prediction of side chains close to helix.\n'
                        'LOOP: Threshold for inclusion of residues for side-chain refinement. Any residues with atoms '
                        'within this distance are included')
    group.addParam('maxEnergyGap', FloatParam, label='Cutoff distance: ',
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
                   help=f'Add terms to the scoring function. You can include multiple TERM keywords to define the '
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
                   help='Click on the wizard to save the specified residues in the list below')
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
                   help='Comma-separated list of residues to add to the refinement list. These should be residues '
                        'that lie outside the distance cutoff')
    group.addParam('selResidue', StringParam, default='',
                   label='Chain of residues: ',
                   condition=f'(stageType == {COMPILE} and residuesDefinition == 1) or (stageType in [{PHELIX}, {PLOOP}]) or (stageType == {TRIM} and trimResidues == 0)',
                   help='Comma-separated list of residues to add to the refinement list. These should be residues '
                        'that lie outside the distance cutoff')

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
      if not pName in msjDic:
        msjDic[pName] = paramDic[pName].default
    return msjDic

  def createMSJDic(self):
    msjDic = {}
    for pName in self.getStageParamsDic(type='Normal').keys():
      if hasattr(self, pName):
        msjDic[pName] = getattr(self, pName).get()
      else:
        print('Something is wrong with parameter ', pName)

    for pName in self.getStageParamsDic(type='Enum').keys():
      if hasattr(self, pName):
        msjDic[pName] = self.getEnumText(pName)
      else:
        print('Something is wrong with parameter ', pName)
    return msjDic

  ############# UTILS
  def buildSimulateStr(self, msjDic):
    '''Checks the values stored in the msjDic and trnaslates them into msjStr.
        If a value is not found in the msjDic, the default is used'''
    msjDic = self.addDefaultForMissing(msjDic)

    glueArg = '[]'
    if msjDic['glueSolute']:
      glueArg = 'solute'

    # NearT and farT must be at least the boundT
    msjDic['nearT'] = max(msjDic['bondedT'], msjDic['nearT'])
    msjDic['farT'] = max(msjDic['bondedT'], msjDic['farT'])
    timeStepArg = TIMESTEP % (msjDic['bondedT'], msjDic['nearT'], msjDic['farT'])

    pressureArg, barostatArg = '', ''
    method = self._thermoDic[msjDic['thermostat']]
    ensemType = msjDic['ensemType']
    if ensemType not in ['NVE', 'NVT', 'Minimization (Brownian)']:
      pressureArg = PRESSURE % (msjDic['pressure'], msjDic['coupleStyle'].lower())
      barostatArg = BAROSTAT % (msjDic['presMDCons'])
      method = self._baroDic[msjDic['barostat']]

    tensionArg, brownianArg = '', ''
    if ensemType == 'Minimization (Brownian)':
      ensemType = 'NVT'
      method = 'Brownie'
      brownianArg = BROWNIAN % (msjDic['deltaMax'])
    elif ensemType == 'NPgT':
      tensionArg = TENSION % msjDic['surfTension']

    restrainArg = ''
    if msjDic['restrains'] != 'None':
      restrainArg = RESTRAINS % (msjDic['restrains'].lower(), msjDic['restrainForce'])

    if not msjDic['annealing']:
      annealArg = 'off'
      tempArg = msjDic['temperature']
    else:
      annealArg = 'on'
      tempArg = self.parseAnnealing(msjDic['annealTemps'])

    msj_str = MSJ_SYSMD_SIM % (annealArg, os.path.abspath(self._getTmpPath()),
                               glueArg, msjDic['simTime'], timeStepArg, tempArg, pressureArg,
                               tensionArg, ensemType, method, msjDic['tempMDCons'], barostatArg, brownianArg,
                               restrainArg, msjDic['velResamp'], msjDic['trajInterval'])
    return msj_str

  def buildMSJ_str(self):
    # todo: write the file neccessary for IFD
    '''Build the .msj (file used by IFD to specify the jobs performed by Schrodinger)
        defining the input parameters'''
    msj_str = MSJ_SYSMD_INIT

    if self.workFlowSteps.get() in ['', None]:
      msjDic = self.createMSJDic()
      msj_str += self.buildSimulateStr(msjDic)
    else:
      workSteps = self.workFlowSteps.get().split('\n')
      if '' in workSteps:
        workSteps.remove('')
      for wStep in workSteps:
        msjDic = eval(wStep)
        msj_str += self.buildSimulateStr(msjDic)

    return msj_str

  def createGUISummary(self):
    with open(self._getExtraPath("summary.txt"), 'w') as f:
      if self.workFlowSteps.get():
        f.write(self.createSummary())
      else:
        f.write(self.createSummary(self.createMSJDic()))

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


  def getJobName(self):
    files = os.listdir(self._getTmpPath())
    for f in files:
      if f.endswith('.msj'):
        return f.replace('.msj', '')

  def getSchJobId(self):
    jobId = None
    jobListFile = os.path.abspath(self._getTmpPath('jobList.txt'))
    if self.getJobName():
      check_call(jobControlProg + ' -list {} | grep {} > {}'.
                 format(self.getJobName(), self.getJobName(), jobListFile), shell=True)
      with open(jobListFile) as f:
        jobId = f.read().split('\n')[0].split()[0]
    return jobId

  def setAborted(self):
    super().setAborted()
    jobId = self.getSchJobId()
    if jobId:
      print('Killing job: {} with jobName {}'.format(jobId, self.getJobName()))
      check_call(jobControlProg + ' -kill {}'.format(jobId), shell=True)

  def findLastCheckPoint(self):
    CP_files = glob.glob(self._getTmpPath('*_checkpoint'))
    if CP_files:
      return natural_sort(CP_files)[-1]
    return CP_files

  def findLastDataTgz(self):
    data_files = glob.glob(self._getTmpPath('*-out.tgz'))
    if data_files:
      return natural_sort(data_files)[-1]
    return data_files

  def getCurrentJobName(self):
    msj_file = glob.glob(self._getTmpPath('simulation*.msj'))
    if msj_file:
      return os.path.basename(msj_file[-1]).replace('.msj', '')
    else:
      return 'simulation_' + str(rd.randint(1000000, 9999999))






