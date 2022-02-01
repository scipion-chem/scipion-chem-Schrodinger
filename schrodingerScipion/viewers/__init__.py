from pyworkflow.gui.browser import FileTreeProvider
from .viewers_data import SchrodingerDataViewer, MaestroFileHandler
from .viewer_protocol_sitemap import ViewerSitemap
from .viewer_grid import BBoxViewer
from .viewer_docking import ProtGlideDockingViewer
from .viewer_simulation import DesmondSimulationViewer

register = FileTreeProvider.registerFileHandler
register(MaestroFileHandler(), '.mae', '.maegz')