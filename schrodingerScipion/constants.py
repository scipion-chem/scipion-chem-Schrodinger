# coding: latin-1
# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano
# *
# * [1] uam, madrid, Spain
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

SCHRODINGER_HOME = 'SCHRODINGER_HOME' #Acceso by Plugin.getHome()

ATTRIBUTES_MAPPING = {'SiteScore': 'score', 'Dscore': 'druggability', 'size': 'nPoints',
                      'volume': 'volume', 'exposure': 'exposure', 'enclosure': 'enclosure',
                      'contact': 'contact', 'phobic': 'hidrophobic', 'philic': 'hidrophilic',
                      'balance': 'balance', 'don/acc': 'don/acc', 'class': 'class',
                      'contactAtoms': 'contactAtoms', 'contactResidues': 'contactResidues'}

#################################### SYSTEM PREPARATION ############

MSJ_SYSPREP = '''task {
  task = "desmond:auto"
}

build_geometry {
  #add counterion
  %s  
  box = {
     shape = %s
     size = %s
     size_type = %s
  }
  override_forcefield = %s
  rezero_system = %s
  #add salt
  %s  
  #add solvent
  %s
}

assign_forcefield {
  forcefield = %s
}

# command example:
# $SCHRODINGER/utilities/multisim -HOST <hostname> -JOBNAME desmond_trial -m desmond_trial.msj \
desmond_trial.mae -o desmond_trial.cms
'''

ADD_COUNTERION = '''add_counterion = {
     ion = %s
     number = %s
  }'''

ADD_SALT = '''salt = {
     concentration = %s
     negative_ion = %s
     positive_ion = %s
  }'''

SOLVENT = 'solvent = %s'

SIZE_SINGLE = '%s'
SIZE_LIST = '[%s %s %s]'
ANGLES = '[%s %s %s %s %s %s]'

################################ SYSTEM RELAXATION ##################

MSJ_SYSRELAX = '''task {
  task = "desmond:auto"
  set_family = {
    desmond = {
      checkpt.write_last_step = no
    }
  }
}

simulate {
  annealing   = off
  dir         = %s
  glue        = %s
  time        = %s
  timestep    = [0.001 0.001 0.003]
  temperature = %s
  #Pressure
  %s
  #Tension
  %s
  ensemble = {
    class  = %s
    method = %s
    thermostat.tau = %s
    #barostat tau
    %s 
    #Browian delta max
    %s
  }
  
  #restrains
  %s   
  
  randomize_velocity.interval = %s
  eneseq.interval             = 0.3

  trajectory {
  center = solute
  first = 0.0
  format = dtr
  interval = %s
  periodicfix = true
  }
}

# command example:
# $SCHRODINGER/utilities/multisim -HOST <hostname> -JOBNAME desmond_trial -m desmond_trial.msj desmond_trial.mae 
#-o desmond_trial.cms 
'''

PRESSURE = '''pressure = [%s %s]'''
TENSION = '''surface_tension = %s'''
BAROSTAT = '''barostat.tau = %s'''

BROWNIAN = '''brownie.delta_max = %s'''

RESTRAINS = '''restrain    = { atom = %s force_constant = %s }'''


