"""The gsnap-indexing command."""
from .base import Base

class Gsnap_indexing(Base):
    """Make GSNAP index"""

    def run(self):
        from .gsnap_indexing_functions import *

        # Unbuffer stdout and redirect stderr into stdout.  This helps observe logged events in realtime.
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
        os.dup2(sys.stdout.fileno(), sys.stderr.fileno())

        # execute the pipeline stage
        make_index(True)
