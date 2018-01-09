"""The non-host-alignment command."""
from .base import Base

class Non_host_alignment(Base):
    """Perform non-host alignment"""

    def run(self):
        from .non_host_alignment_functions import *

        # Unbuffer stdout and redirect stderr into stdout.  This helps observe logged events in realtime.
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
        os.dup2(sys.stdout.fileno(), sys.stderr.fileno())

        # execute the pipeline stage
        run_stage2(True)
