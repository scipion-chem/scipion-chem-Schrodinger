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

import os

from pyworkflow.protocol.params import PointerParam, StringParam,\
  EnumParam, BooleanParam, FloatParam, IntParam, LEVEL_ADVANCED
from pwem.protocols import EMProtocol
from schrodingerScipion import Plugin as schrodinger_plugin
from schrodingerScipion.constants import *
from schrodingerScipion.objects import SchrodingerAtomStruct

multisimProg = schrodinger_plugin.getHome('utilities/multisim')

STRUCTURE, LIGAND = 0, 1

class ProtSchrodingerDesmondSysRelax(EMProtocol):
    """Calls Desmond molecular dynamics for the preparation of the system via solvatation, the addition of ions
    and a force field"""
    _label = 'system relaxation (desmond)'
    _program = ""

    _relaxTypes = ['NVE', 'NVT', 'NPT', 'NPAT', 'NPgT', 'Minimization (Brownian)']
    _thermoDic = {'Noose-Hover': 'NH', 'Langevin': 'Langevin', 'None': 'None'}
    _thermostats = ['Noose-Hover', 'Langevin', 'Brownie', 'None']

    _baroDic = {'Martyna-Tobias-Klein': 'MTK', 'Langevin': 'Langevin', 'None': 'None'}
    _barostats = ['Martyna-Tobias-Klein', 'Langevin', 'None']
    _coupleStyle = ['Isotropic', 'Semi-isotropic', 'Anisotropic', 'Constant area']
    _restrainTypes = ['None', 'Ligand', 'Protein', 'Solute_heavy_atom', 'Solute']

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputStruct', PointerParam, pointerClass='SchrodingerAtomStruct',
                      label='Input system to be prepared for MD:', allowsNull=False,
                      help='Atomic structure to be relaxed')
        group = form.addGroup('Simulation time')
        group.addParam('simTime', FloatParam, default=100.0,
                       label='Relaxation time (ps):',
                       help='Time of the simulation step (ns)')
        group.addParam('velResamp', FloatParam, default=1.0,
                       label='Velocity resampling (ps):', expertLevel=LEVEL_ADVANCED,
                       help='Velocities of the particles are resampled every x ps')

        group = form.addGroup('Trajectory')
        group.addParam('glueSolute', BooleanParam, default=True,
                      label='Glue close solute molecules: ', expertLevel=LEVEL_ADVANCED,
                      help='Glue close solute molecules together. Does not affect energy, only for visualization'
                           ' of the trajectory')
        group.addParam('trajInterval', FloatParam, default=5.0,
                       label='Interval time (ps):',
                       help='Time between each frame recorded in the simulation (ns)')


        form.addSection('Relaxation')
        group = form.addGroup('Ensemble')
        group.addParam('ensemType', EnumParam, default=0,
                       label='Relaxation type: ', choices=self._relaxTypes,
                       help='Relax type of the simulation')

        line = group.addLine('Temperature settings: ', condition='ensemType!=0',
                             help='Temperature during the simulation (K)\nThermostat type\n'
                                  'Relaxation time constant for thermostat (ps)')
        line.addParam('temperature', FloatParam, default=300.0,
                       label='Temperature: ')
        line.addParam('thermostat', EnumParam, default=0, condition='ensemType!=5',
                       label='Thermostat type: ', choices=self._thermostats, expertLevel=LEVEL_ADVANCED)
        line.addParam('deltaMax', FloatParam, default=0.1, condition='ensemType==5',
                      label='Max displacement: ', expertLevel=LEVEL_ADVANCED)
        line.addParam('tempRelaxCons', FloatParam, default=1.0,
                       label='Temperature relax constant: ', expertLevel=LEVEL_ADVANCED)

        line = group.addLine('Pressure settings: ', condition='ensemType not in [0, 1, 5]',
                             help='Pressure during the simulation (bar)\nBarostat type\n'
                                  'Relaxation time constant for barostat (ps)')
        line.addParam('pressure', FloatParam, default=1.01325,
                       label='   Pressure:   ')
        line.addParam('barostat', EnumParam, default=0,
                       label='  Barostat type:   ', choices=self._barostats, expertLevel=LEVEL_ADVANCED)
        line.addParam('presRelaxCons', FloatParam, default=2.0,
                       label='   Pressure relax constant:   ', expertLevel=LEVEL_ADVANCED)
        line = group.addLine('Pressure coupling: ', condition='ensemType not in [0, 1, 5]',
                             help='Pressure coupling style', expertLevel=LEVEL_ADVANCED)
        line.addParam('coupleStyle', EnumParam, default=0,
                      label='  Coupling style:   ', choices=self._coupleStyle)

        line = group.addLine('Tension settings: ', condition='ensemType==4',
                             help='Surface tension during the simulation (bar·Å)')
        line.addParam('surfTension', FloatParam, default=0.0,
                      label='Surface tension: ')

        group = form.addGroup('Restrains')
        group.addParam('restrains', EnumParam, default=0,
                       label='Restrains: ', choices=self._restrainTypes,
                       help='Restrain movement of specific groups of atoms')
        group.addParam('restrainForce', FloatParam, default=50.0,
                       label='Restrain force constant: ', condition='restrains!=0',
                       help='Restrain force applied to the selection (kcal/mol/Å2)')


    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('relaxationStep')

    def relaxationStep(self):
        maeFile = self.inputStruct.get().getFileName()
        sysName = maeFile.split('/')[-1].split('.')[0]

        msjFile = self._getExtraPath('{}.msj'.format(sysName))
        msjStr = self.buildMSJ_str()
        with open(msjFile, 'w') as f:
            f.write(msjStr)

        cmsFile = sysName+'-out.cms'

        args = ' -m {} {} -o {} -WAIT -JOBNAME {}'.format(msjFile.split('/')[-1], os.path.abspath(maeFile),
                                                          cmsFile, sysName)
        self.runJob(multisimProg, args, cwd=self._getExtraPath())

        if os.path.exists(cmsFile):
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
        glueArg = '[]'
        if self.glueSolute.get():
            glueArg = 'solute'

        pressureArg, barostatArg = '', ''
        method = self._thermoDic[self.getEnumText('thermostat')]
        ensemType = self.getEnumText('ensemType')
        if self.ensemType.get() not in [0, 1, 5]:
            pressureArg = PRESSURE % (self.pressure.get(), self.getEnumText('coupleStyle').lower())
            barostatArg = BAROSTAT % (self.presRelaxCons.get())
            method = self._baroDic[self.getEnumText('barostat')]

        tensionArg = ''
        if self.ensemType.get() == 5:
            ensemType = 'NVT'
            method = 'Brownie'
        elif self.ensemType.get() == 4:
            tensionArg = TENSION % self.surfTension.get()

        brownianArg = ''
        if self.ensemType.get() == 5:
            brownianArg = BROWNIAN % (self.deltaMax.get())

        restrainArg = ''
        if self.restrains.get() != 0:
            restrainArg = RESTRAINS % (self.getEnumText('restrains').lower(), self.restrainForce.get())

        msj_str = MSJ_SYSRELAX % (os.path.abspath(self._getExtraPath()), glueArg, self.simTime.get(),
                                  self.temperature.get(), pressureArg, tensionArg, ensemType, method,
                                  self.tempRelaxCons.get(), barostatArg, brownianArg, restrainArg,
                                  self.velResamp.get(), self.trajInterval.get())
        return msj_str


