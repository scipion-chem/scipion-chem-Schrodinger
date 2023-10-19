# **************************************************************************
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
# General imports
import os

# Scipion em imports
from pwem.viewers.views import ObjectView
import pyworkflow.viewer as pwviewer
from pyworkflow.gui.browser import FileHandler
import pyworkflow.utils as pwutils

# Plugin imports
from .. import Plugin
from ..objects import SchrodingerAtomStruct, SchrodingerBindingSites, SchrodingerSystem

class SchrodingerDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        SchrodingerAtomStruct,
        SchrodingerBindingSites,
        SchrodingerSystem
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        # For now handle both types of SetOfTiltSeries together
        if issubclass(cls, SchrodingerAtomStruct) or \
           issubclass(cls, SchrodingerBindingSites) or \
           issubclass(cls, SchrodingerSystem):
            cwd = '/'.join(obj.getFileName().split('/')[:-1])
            views.append(MaestroView(os.path.abspath(obj.getFileName()), cwd=cwd))
        return views

class MaestroView(pwviewer.CommandView):
    """ View for calling an external command. """
    def __init__(self, inputFile, **kwargs):
        pwviewer.CommandView.__init__(self, "%s %s &"%(Plugin.getHome('maestro'),inputFile),
                                      env=Plugin.getEnviron(), **kwargs)

class MaestroFileHandler(FileHandler):
    def getFileActions(self, objFile):
        fn = objFile.getPath()
        return [('Open with Maestro', lambda: MaestroView(fn).show(),
                 pwutils.Icon.ACTION_VISUALIZE)]

    def getFileIcon(self, objFile):
        return 'file_vol.gif'

