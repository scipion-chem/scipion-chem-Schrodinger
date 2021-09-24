from pyworkflow.gui.browser import FileTreeProvider
from .viewers_data import SchrodingerDataViewer, MaestroFileHandler
from .viewer_glide_docking import ProtSchrodingerGlideDockingViewer
from .viewer_protocol_sitemap import ViewerSitemap
from .viewer_grid import BBoxViewer

register = FileTreeProvider.registerFileHandler
register(MaestroFileHandler(), '.mae', '.maegz')