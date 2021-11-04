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

import os, shutil, threading, subprocess

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, EnumParam, BooleanParam, FloatParam, IntParam, STEPS_PARALLEL
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.utils.path import makePath
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import relabelAtomsMol2
from schrodingerScipion import Plugin as schrodinger_plugin
from schrodingerScipion.constants import *
from schrodingerScipion.utils.utils import putMol2Title, sortDockingResults
from schrodingerScipion.objects import SchrodingerPoses, SchrodingerAtomStruct

multisimProg = schrodinger_plugin.getHome('utilities/multisim')
glideProg = schrodinger_plugin.getHome('glide')
structConvertProg = schrodinger_plugin.getHome('utilities/structconvert')
structCatProg = schrodinger_plugin.getHome('utilities/structcat')
propListerProg = schrodinger_plugin.getHome('utilities/proplister')
maeSubsetProg = schrodinger_plugin.getHome('utilities/maesubset')

STRUCTURE, LIGAND = 0, 1


class ProtSchrodingerDesmondSysPrep(EMProtocol):
    """Calls Desmond molecular dynamics for the preparation of the system via solvatation, the addition of ions
    and a force field"""
    _label = 'system preparation (desmond)'
    _program = ""

    _solventTypes = ['SPC', 'TIP3P', 'TIP4P', 'TIP4PEW', 'TIP4PD', 'TIP5P',
                     'SPCE', 'TIP4P2005', 'DMSO', 'METHANOL', 'OCTANOL']
    _boundaryShapes = ['Cubic', 'Orthorhombic', 'Triclinic',
                       'Truncated octahedron', 'Rhombic dodecahedron xy-square', 'Rhombic dodecahedron xy-hexagon']
    _cations = ['Na+', 'Li+', 'K+', 'Rb+', 'Cs+', 'Mg2+', 'Ca2+', 'Zn2+', 'Fe2+', 'Fe3+']
    _anions = ['F-', 'Cl-', 'Br-', 'I-']

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputFrom', EnumParam, default=STRUCTURE,
                      label='Input from: ', choices=['AtomStructure', 'SetOfSmallMolecules'],
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

        form.addParam('rezero', BooleanParam, default=True,
                       label='Reset origin to center of mass: ',
                       help='Set system origin as the solute center of mass')

        group = form.addGroup('Boundary box')
        group.addParam('boundaryShape', EnumParam,
                      choices=self._boundaryShapes, default=1,
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
                      help='Minimize volume of the resulting box')

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
        group.addParam('placeIons', EnumParam, default=0,
                       label='Add ions: ', choices=['None', 'Neutralize', 'Add number'],
                       help='Whether to add ions to the system')
        line = group.addLine('Ion type:', condition='placeIons!=0',
                             help='Type of the ions to neutralize charges (Depending on the system charge)')
        line.addParam('cationType', EnumParam, condition='placeIons!=0',
                       label='Cation to add: ', choices=self._cations,
                       help='Whether to add cations to the system')
        line.addParam('anionType', EnumParam, condition='placeIons!=0',
                       label='Anion to add: ', choices=self._anions,
                       help='Whether to add anions to the system')
        group.addParam('ionNum', IntParam, condition='placeIons==2',
                       label='Number of ions to add: ',
                       help='Number of ions to add, cations if the system charge is negative and viceversa')

        group = form.addGroup('Salt')
        group.addParam('addSalt', BooleanParam, default=False,
                       label='Add a salt into the system: ',
                       help='Add a salt into the system')
        group.addParam('saltConc', FloatParam, condition='addSalt',
                       default=0.15,
                       label='Salt concentration: ',
                       help='Salt concentration')
        line = group.addLine('Salt type:', condition='addSalt',
                             help='Type of the ions to neutralize charges (Depending on the system charge)')
        line.addParam('cationTypeSalt', EnumParam, condition='addSalt',
                      label='Salt cation: ', choices=self._cations,
                      help='Cation type of the salt')
        line.addParam('anionTypeSalt', EnumParam, condition='addSalt',
                       label='Salt anion: ', choices=self._anions,
                       help='Anion type of the salt')

        group = form.addGroup('Force field')
        group.addParam('ffType', EnumParam,
                      label='Force field: ', choices=['S-OPLS', 'OPLS_2005'], default=0,
                      help='Force field to use')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('solutePreparationStep')
        self._insertFunctionStep('systemPreparationStep')

    def solutePreparationStep(self):
        if self.inputFrom.get() == STRUCTURE:
            if type(self.inputStruct.get()) == SchrodingerAtomStruct:
                self.soluteFile = self.inputStruct.get().getFileName()
            else:
                pdbFile = self.inputStruct.get().getFileName()
                structName = os.path.splitext(os.path.basename(pdbFile))[0]
                self.soluteFile = self._getExtraPath(structName + '.mae')
                self.runJob(structConvertProg, '{} {}'.format(pdbFile, self.soluteFile))

        elif self.inputFrom.get() == LIGAND:
            mol = self.getSpecifiedMol()
            molMaeFile = self._getExtraPath(mol.getUniqueName() + '.maegz')
            self.runJob(structConvertProg, '{} {}'.format(mol.getPoseFile(), molMaeFile))

            if hasattr(mol, 'structFile'):
                targetMaeFile = mol.structFile
            else:
                targetFile = self.inputSetOfMols.get().getProteinFile()
                targetName = os.path.splitext(os.path.basename(targetFile))[0]
                targetMaeFile = self._getExtraPath(targetName + '.maegz')
                self.runJob(structConvertProg, '{} {}'.format(targetFile, targetMaeFile))

            self.soluteFile = self._getExtraPath('complexSolute.maegz')
            self.runJob('zcat', '{} {} > {}'.format(molMaeFile, targetMaeFile, self.soluteFile))


    def systemPreparationStep(self):
        maeFile = self.soluteFile
        sysName = maeFile .split('/')[-1].split('.')[0]

        msjFile = self._getExtraPath('{}.msj'.format(sysName))
        msjStr = self.buildMSJ_str()
        with open(msjFile, 'w') as f:
            f.write(msjStr)

        cmsFile = sysName+'-out.cms'

        args = ' -m {} {} -o {} -WAIT -JOBNAME {}'.format(msjFile.split('/')[-1], os.path.abspath(maeFile),
                                                          cmsFile, sysName)
        self.runJob(multisimProg, args, cwd=self._getExtraPath())

        cmsStruct = SchrodingerAtomStruct()
        cmsStruct.setFileName(self._getExtraPath(cmsFile))
        self._defineOutputs(outputSystem=cmsStruct)


    def _validate(self):
        errors = []
        return errors

    ############# UTILS

    def buildMSJ_str(self):
        '''Build the .msj (file used by multisim to specify the jobs performed by Schrodinger)
        defining the input parameters'''
        addIonsArg = ''
        if self.placeIons.get() != 0:
            #todo: check if cation or anion
            if self.placeIons.get() == 1:
                number = 'neutralize_system'
            else:
                number = self.ionNum.get()
            addIonsArg = ADD_COUNTERION % (self.getEnumText('cationType')[:-1], number)

        boxArgs = [self.getEnumText('boundaryShape').lower()]
        if self.boundaryShape.get() in [1, 2]:
            boxArgs += [self.distA.get(), self.distB.get(), self.distC.get()]
        else:
            boxArgs += [self.distA.get(), '', '']

        if self.boundaryShape.get() == 2:
            boxArgs += [ANGLES % (self.angleA.get(), self.angleB.get(), self.angleC.get())]
        else:
            boxArgs += ['']

        boxArgs += [self.getEnumText('boxMethod').lower()]

        saltArg = ''
        if self.addSalt:
            saltArg = ADD_SALT % (self.saltConc.get(), self.getEnumText('anionTypeSalt')[:-1],
                                  self.getEnumText('cationTypeSalt')[:-1])

        solventArg = ''
        if self.solvate:
            solventArg = SOLVENT % self.getEnumText('solventType')

        msj_str = MSJ_SYSPREP % (addIonsArg, *boxArgs, self.getEnumText('ffType'), str(self.rezero.get()).lower(),
                                 saltArg, solventArg, self.getEnumText('ffType'))
        return msj_str


    def getSpecifiedMol(self):
        myMol = None
        for mol in self.inputSetOfMols.get():
          if mol.getUniqueName() == self.inputLigand.get():
            myMol = mol.clone()
            break
        if myMol == None:
            print('The input ligand is not found')
            return None
        else:
            return myMol


