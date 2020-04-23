from pyworkflow.gui.browser import FileTreeProvider
from .viewers_data import SchrodingerDataViewer, MaestroFileHandler
from .viewer_glide_docking import ProtSchrodingerGlideDockingViewer

register = FileTreeProvider.registerFileHandler
register(MaestroFileHandler(), '.mae', '.maegz')