# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

"""
This package contains protocols interfacing Schrodinger's Maestro
"""

import os
import pwem
import pyworkflow.utils as pwutils
from .bibtex import _bibtexStr

_logo = 'schrodinger.png'

class Plugin(pwem.Plugin):
    _homeVar = 'SCHRODINGER_HOME'

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar('SCHRODINGER_HOME', 'schrodinger2019-4')

    @classmethod
    def defineBinaries(cls, env):
        pass

    @classmethod
    def getEnviron(cls, schrodingerFirst=True):
        """ Create the needed environment for Schrodinger programs. """
        environ = pwutils.Environ(os.environ)
        pos = pwutils.Environ.BEGIN if schrodingerFirst else pwutils.Environ.END
        environ.update({
            'PATH': cls.getHome('')
        }, position=pos)

        return environ
