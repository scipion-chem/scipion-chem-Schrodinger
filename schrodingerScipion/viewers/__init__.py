from pyworkflow.gui.browser import FileTreeProvider
from .viewers_data import SchrodingerDataViewer, MaestroFileHandler

register = FileTreeProvider.registerFileHandler
register(MaestroFileHandler(), '.mae', '.maegz')