"""The gsnap-indexing command."""
from .base import Base


class Gsnap_indexing(Base):
    """Make GSNAP index"""

    def run(self):
        from .gsnap_indexing_functions import *

        set_up_stdout()
        make_index(self.version)
