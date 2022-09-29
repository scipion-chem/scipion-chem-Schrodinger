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

import os, time, subprocess
import random as rd
from subprocess import check_call

from pyworkflow.protocol.params import PointerParam, StringParam,\
  EnumParam, BooleanParam, FloatParam, IntParam, LEVEL_ADVANCED
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toPdb

from pwchem.utils import pdbqt2other, convertToSdf

from schrodingerScipion import Plugin as schrodingerPlugin
from schrodingerScipion.constants import *
from schrodingerScipion.objects import SchrodingerAtomStruct, SchrodingerSystem
from schrodingerScipion.utils import getChargeFromMAE

multisimProg = schrodingerPlugin.getHome('utilities/multisim')
jobControlProg = schrodingerPlugin.getHome('jobcontrol')
structConvertProg = schrodingerPlugin.getHome('utilities/structconvert')
progLigPrep = schrodingerPlugin.getHome('ligprep')

STRUCTURE, LIGAND = 0, 1

class ProtSchrodingerDesmondSysPrep(EMProtocol):
    """Calls Desmond molecular dynamics for the preparation of the system via solvatation, the addition of ions
    and a force field"""
    _label = 'system preparation (desmond)'
    _program = ""

    _solventTypes = ['SPC', 'TIP3P', 'TIP4P', 'TIP4PEW', 'TIP4PD', 'TIP5P',
                     'DMSO', 'METHANOL', 'OCTANOL']
    _boundaryShapes = {'Cubic': 'cubic', 'Orthorhombic': 'orthorhombic', 'Triclinic': 'triclinic',
                       'Truncated octahedron': 'truncated_octahedron',
                       'Rhombic dodecahedron xy-square': 'dodecahedron_square',
                       'Rhombic dodecahedron xy-hexagon': 'dodecahedron_hexagon'}
    _saltCations = ['Na+', 'Li+', 'K+', 'Rb+', 'Cs+', 'Mg2+', 'Ca2+', 'Zn2+', 'Fe2+', 'Fe3+']
    _cations = ['Na+', 'Li+', 'K+', 'Rb+', 'Cs+']
    _anions = ['F-', 'Cl-', 'Br-', 'I-']

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputFrom', EnumParam, default=STRUCTURE,
                      label='Input from: ', choices=['AtomStruct', 'SetOfSmallMolecules'],
                      help='Type of input you want to use')
        form.addParam('inputStruct', PointerParam, pointerClass='SchrodingerAtomStruct, AtomStruct',
                      label='Input structure to be prepared for MD:', allowsNull=False, condition='inputFrom==0',
                      help='Atomic structure to be prepared for MD by solvation, ions addition etc')
        form.addParam('inputSetOfMols', PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Input set of molecules:', allowsNull=False, condition='inputFrom==1',
                      help='Input set of docked molecules. One of them will be prepared together with its target')
        form.addParam('inputLigand', StringParam, condition='inputFrom==1',
                      label='Ligand to prepare: ',
                      help='Specific ligand to prepare in the system')

        form.addParam('prepareTarget', BooleanParam, default=True,
                      label='Prepare target: ', expertLevel=LEVEL_ADVANCED,
                      help='Prepare target with Schrodinger PrepWizard to ensure correct structure.'
                           'It may subtly variate the atom positions')
        form.addParam('prepareLigand', BooleanParam, default=True, condition='inputFrom=={}'.format(LIGAND),
                      label='Prepare ligand: ', expertLevel=LEVEL_ADVANCED,
                      help='Prepare target with Schrodinger LigPrep to ensure correct structure.'
                           'It may subtly variate the ligand atom positions')

        group = form.addGroup('Boundary box')
        group.addParam('boundaryShape', EnumParam,
                      choices=list(self._boundaryShapes.keys()), default=1,
                      label='Solvent type: ',
                      help='Different water and other chemics models to use as solvent')
        group.addParam('boxMethod', EnumParam,
                       choices=['Absolute', 'Buffer'], default=1,
                       label='Box size method: ',
                       help='Whether to use absolute size values or minimum distances to the '
                            'structure to build the box')
        line = group.addLine('Box size:',
                            help='Distances of the bounding box')
        line.addParam('distA', FloatParam,
                       default=10.0, label='A: ')
        line.addParam('distB', FloatParam, condition='boundaryShape in [1, 2]',
                       default=10.0, label='B: ')
        line.addParam('distC', FloatParam, condition='boundaryShape in [1, 2]',
                       default=10.0, label='C: ')
        line = group.addLine('Angles:', condition='boundaryShape == 2',
                             help='Angles of the bounding box')
        line.addParam('angleA', FloatParam, condition='boundaryShape == 2',
                       default=90.0, label='A: ')
        line.addParam('angleB', FloatParam, condition='boundaryShape == 2',
                       default=90.0, label='B: ')
        line.addParam('angleC', FloatParam, condition='boundaryShape == 2',
                       default=90.0, label='C: ')
        group.addParam('minimize', BooleanParam, default=False,
                      label='Minimize volume: ',
                      help='Minimize volume of the resulting box by rotating the solute to fit better in the box')

        group = form.addGroup('Solvation model')
        group.addParam('solvate', BooleanParam, default=True,
                      label='Solvate the atomic structure: ',
                      help='Introduce the structure into a box with a solvent')
        group.addParam('solventType', EnumParam, condition='solvate',
                      choices=self._solventTypes, default=0,
                      label='Solvent type: ',
                      help='Different water and other chemics models to use as solvent')

        form.addSection('Charges')
        group = form.addGroup('Ions')
        group.addParam('placeIons', EnumParam, default=1,
                       label='Add ions: ', choices=['None', 'Neutralize', 'Add number'],
                       help='Whether to add ions to the system')
        line = group.addLine('Solute charge:', condition='placeIons!=0', expertLevel=LEVEL_ADVANCED)
        line.addParam('solCharge', IntParam, default=0, condition='placeIons!=0', readOnly=True,
                      help='Check charge of the solute before the addition of ions')
        line = group.addLine('Cation type:', condition='placeIons!=0',
                             help='Type of the cations to add into the system. '
                                  '(If neutralize, only added when solute has negative charge)')
        line.addParam('cationType', EnumParam, condition='placeIons!=0',
                       label='Cation to add: ', choices=self._cations, default=0,
                       help='Which cations to add in the system')
        line.addParam('cationNum', IntParam, condition='placeIons==2',
                      label='Number of cations to add: ',
                      help='Number of cations to add')

        line = group.addLine('Anion type:', condition='placeIons!=0',
                             help='Type of the anions to add into the system. '
                                  '(If neutralize, only added when solute has positive charge)')
        line.addParam('anionType', EnumParam, condition='placeIons!=0',
                      label='Anion to add: ', choices=self._anions, default=1,
                      help='Which anions to add in the system')
        line.addParam('anionNum', IntParam, condition='placeIons==2',
                      label='Number of anions to add: ',
                      help='Number of anions to add')

        group = form.addGroup('Salt')
        group.addParam('addSalt', BooleanParam, default=False,
                       label='Add a salt into the system: ',
                       help='Add a salt into the system')
        group.addParam('saltConc', FloatParam, condition='addSalt',
                       default=0.15,
                       label='Salt concentration (M): ',
                       help='Salt concentration')
        line = group.addLine('Salt type:', condition='addSalt',
                             help='Type of the ions to neutralize charges (Depending on the system charge)')
        line.addParam('cationTypeSalt', EnumParam, condition='addSalt',
                      label='Salt cation: ', choices=self._saltCations, default=0,
                      help='Cation type of the salt')
        line.addParam('anionTypeSalt', EnumParam, condition='addSalt',
                       label='Salt anion: ', choices=self._anions, default=1,
                       help='Anion type of the salt')

        group = form.addGroup('Force field')
        group.addParam('ffType', EnumParam,
                      label='Force field: ', choices=['S-OPLS', 'OPLS_2005'], default=0,
                      help='Force field to use')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('solutePreparationStep')
        self._insertFunctionStep('systemPreparationStep')

    def solutePreparationStep(self):
        if self.inputFrom.get() == STRUCTURE:
            if type(self.inputStruct.get()) == SchrodingerAtomStruct:
                self.soluteFile = self.inputStruct.get().getFileName()
            else:
                pdbFile = self.getPdbFile()
                structName = os.path.splitext(os.path.basename(pdbFile))[0]
                self.soluteFile = os.path.abspath(self._getExtraPath(structName + '.mae'))
                if not os.path.exists(self.soluteFile):
                  if self.prepareTarget.get():
                      self.soluteFile = self.prepareTargetFile(pdbFile, self.soluteFile)
                  else:
                      self.runJob(structConvertProg, '{} {}'.format(pdbFile, self.soluteFile))

        elif self.inputFrom.get() == LIGAND:
            self.soluteFile = self._getExtraPath('complexSolute.mae')
            if not os.path.exists(self.soluteFile):
                mol = self.getSpecifiedMol()
                molFile = mol.getPoseFile()
                if molFile.endswith('.pdbqt'):
                    molFile = convertToSdf(self, molFile)

                if self.prepareLigand.get():
                    molMaeFile = self.prepareLigandFile(molFile)
                else:
                    molMaeFile = self._getExtraPath(mol.getUniqueName() + '.maegz')
                    self.runJob(structConvertProg, '{} {}'.format(molFile, molMaeFile))

                if hasattr(mol, 'structFile'):
                    targetMaeFile = mol.structFile
                else:
                    pdbFile = self.getPdbFile()
                    targetName = os.path.splitext(os.path.basename(pdbFile))[0]
                    targetMaeFile = os.path.abspath(self._getExtraPath(targetName + '.maegz'))
                    if self.prepareTarget.get():
                        targetMaeFile = self.prepareTargetFile(pdbFile, targetMaeFile)
                    else:
                        self.runJob(structConvertProg, '{} {}'.format(pdbFile, targetMaeFile))

                self.runJob('zcat', '{} {} > {}'.format(molMaeFile, targetMaeFile, self.soluteFile))

    def systemPreparationStep(self):
        maeFile = self.soluteFile
        sysName = maeFile.split('/')[-1].split('.')[0]
        jobName = sysName + '_' + str(rd.randint(1000000, 9999999))

        msjFile = self._getExtraPath('{}.msj'.format(sysName))
        msjStr = self.buildMSJ_str(maeFile)
        with open(msjFile, 'w') as f:
            f.write(msjStr)

        cmsName = sysName+'-out.cms'
        args = ' -m {} {} -o {} -WAIT -JOBNAME {}'.format(msjFile.split('/')[-1], os.path.abspath(maeFile),
                                                          cmsName, jobName)
        self.runJob(multisimProg, args, cwd=self._getExtraPath())

        cmsStruct = SchrodingerSystem()
        if os.path.exists(self._getExtraPath(cmsName)):
            cmsFile = os.path.abspath(self._getPath(cmsName))
            os.rename(self._getExtraPath(cmsName), cmsFile)

            cmsStruct.setFileName(cmsFile)
            self._defineOutputs(outputSystem=cmsStruct)
        else:
            print('Schrodinger system preparation failed')
            open(self._getExtraPath(cmsName))


    def _validate(self):
        errors = []
        return errors

    ############# UTILS

    def buildMSJ_str(self, maeFile):
        '''Build the .msj (file used by multisim to specify the jobs performed by Schrodinger)
        defining the input parameters'''
        addIonsArg = ''
        if self.placeIons.get() != 0:
            solCharge = getChargeFromMAE(maeFile)
            if solCharge < 0:
                number = 'neutralize_system' if self.placeIons.get() == 1 else self.cationNum.get()
                addIonsArg = ADD_COUNTERION % (self.getEnumText('cationType')[:-1], number)
            elif solCharge > 0:
                number = 'neutralize_system' if self.placeIons.get() == 1 else self.anionNum.get()
                addIonsArg = ADD_COUNTERION % (self.getEnumText('anionType')[:-1], number)


        boxArgs = [self._boundaryShapes[self.getEnumText('boundaryShape')]]
        if self.boundaryShape.get() == 1:
            boxArgs += [SIZE_LIST % (self.distA.get(), self.distB.get(), self.distC.get())]

        elif self.boundaryShape.get() == 2:
            boxArgs += [ANGLES % (self.distA.get(), self.distB.get(), self.distC.get(),
                                  self.angleA.get(), self.angleB.get(), self.angleC.get())]
        else:
            boxArgs += [SIZE_SINGLE % self.distA.get()]

        boxArgs += [self.getEnumText('boxMethod').lower()]

        saltArg = ''
        if self.addSalt:
            saltArg = ADD_SALT % (self.saltConc.get(), self.getEnumText('anionTypeSalt')[:-1],
                                  self.getEnumText('cationTypeSalt')[:-1])

        solventArg = ''
        if self.solvate:
            solventArg = SOLVENT % self.getEnumText('solventType')

        msj_str = MSJ_SYSPREP % (addIonsArg, *boxArgs, self.getEnumText('ffType'), "true",
                                 saltArg, solventArg, self.getEnumText('ffType'))
        return msj_str


    def getSpecifiedMol(self):
        myMol = None
        for mol in self.inputSetOfMols.get():
          if mol.__str__() == self.inputLigand.get():
            myMol = mol.clone()
            break
        if myMol == None:
            print('The input ligand is not found')
            return None
        else:
            return myMol

    def getJobName(self):
        files = os.listdir(self._getExtraPath())
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

    def getPdbFile(self):
      if self.inputFrom.get() == STRUCTURE:
          proteinFile = self.inputStruct.get().getFileName()
      elif self.inputFrom.get() == LIGAND:
          proteinFile = self.inputSetOfMols.get().getProteinFile()
      inName, inExt = os.path.splitext(os.path.basename(proteinFile))

      if inExt == '.pdb':
          return os.path.abspath(proteinFile)
      else:
        pdbFile = os.path.abspath(os.path.join(self._getExtraPath(inName + '.pdb')))
        if inExt == '.pdbqt':
            pdbqt2other(self, proteinFile, pdbFile)
        else:
            toPdb(proteinFile, pdbFile)
        return os.path.abspath(pdbFile)

    def prepareLigandFile(self, sdfFile, maeFile=None):
        # Manage files from autodock: 1) Convert to readable by schro (SDF). 2) correct preparation.
        if not sdfFile.endswith('.sdf'):
            sdfFile = convertToSdf(self, sdfFile)

        if not maeFile:
            baseName = os.path.splitext(os.path.basename(sdfFile))[0]
            maeFile = os.path.abspath(os.path.join(self._getExtraPath(baseName + '.maegz')))

        args = " -i 1 -nt -ns -a -isd {} -omae {}".format(sdfFile, maeFile)
        subprocess.check_call([progLigPrep, *args.split()])
        while not os.path.exists(maeFile):
            time.sleep(0.2)
        return maeFile

    def prepareTargetFile(self, inFile, outFile):
        prog = schrodingerPlugin.getHome('utilities/prepwizard')
        args = '-WAIT -noprotassign -noimpref -noepik '
        args += '%s %s' % (inFile, outFile)
        self.runJob(prog, args, cwd=self._getPath())
        return outFile