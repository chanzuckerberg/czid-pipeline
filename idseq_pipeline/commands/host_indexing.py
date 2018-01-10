"""The host-indexing command."""
from .base import Base

class Host_indexing(Base):
    """Perform indexing"""

    def run(self):
        from .host_indexing_functions import *

        # Unbuffer stdout and redirect stderr into stdout.  This helps observe logged events in realtime.
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
        os.dup2(sys.stdout.fileno(), sys.stderr.fileno())

        # make indexes
        make_indexes(True)
