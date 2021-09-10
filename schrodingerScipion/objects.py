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
from pwchem.objects import ProteinPocket
from .utils.utils import parseLogProperties
from .constants import ATTRIBUTES_MAPPING as AM

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

class SitemapPocket(ProteinPocket):
    """ Represent a pocket file from Sitemap"""
    def __init__(self, filename=None, proteinFile=None, logFile=None, **kwargs):
        if filename != None:
            self.pocketId = int(filename.split('-')[1].split('.')[0])
            if logFile != None:
                self.properties = parseLogProperties(logFile, self.pocketId)
                self.properties['class'] = 'SiteMap'
                kwargs.update(self.getKwargs(self.properties, AM))

        super().__init__(filename, proteinFile, **kwargs)
        if hasattr(self, 'pocketId'):
            self.setObjId(self.pocketId)

        #Build contact atoms
        if proteinFile != None:
            cAtoms = self.buildContactAtoms(calculate=True)
            self.setContactAtoms(self.encodeIds(self.getAtomsIds(cAtoms)))
            cResidues = self.getResiduesFromAtoms(cAtoms)
            self.setContactResidues(self.encodeIds(self.getResiduesIds(cResidues)))

    def __str__(self):
        s = 'SiteMap pocket {}\nFile: {}'.format(self.getObjId(), self.getFileName())
        return s
