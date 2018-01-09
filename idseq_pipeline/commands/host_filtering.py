"""The host-filtering command."""
from .base import Base

class Host_filtering(Base):
    """Perform host filtering"""

    def run(self):
        from .host_filtering_functions import *

        # Unbuffer stdout and redirect stderr into stdout.  This helps observe logged events in realtime.
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
        os.dup2(sys.stdout.fileno(), sys.stderr.fileno())

        # execute the pipeline stage
        run_stage1(True)
