"""The lineage command."""
from .base import Base


class Lineages(Base):
    """Make lineage files"""

    def run(self):
        from .lineages_functions import *

        unbuffer_stdout()
        make_lineages(self.version)
