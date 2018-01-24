"""The rapsearch-indexing command."""
from .base import Base

class Rapsearch_indexing(Base):
    """Make RAPSearch2 index"""

    def run(self):
        from .rapsearch_indexing_functions import *

        unbuffer_stdout()
        make_index(self.version)
