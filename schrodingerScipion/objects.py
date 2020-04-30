# -*- coding: utf-8 -*-
#  **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
import pwem.objects.data as data

class SchrodingerAtomStruct(data.EMFile):
    """An AtomStruct in the file format of Maestro"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)

    def getExtension(self):
        return os.path.splitext(self.getFileName())[1]

class SchrodingerGrid(data.EMFile):
    """A search grid in the file format of Maestro"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)

class SchrodingerBindingSites(data.EMFile):
    """A set of binding sites in the file format of Maestro"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)

class SchrodingerPoses(data.EMFile):
    """A set of poses and a structure in the file format of Maestro"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)
