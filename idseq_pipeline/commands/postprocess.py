"""The post-processing command."""
from .base import Base

class Postprocess(Base):
    """Perform post-processing"""

    def run(self):
        from .postprocess_functions import *

        # Unbuffer stdout and redirect stderr into stdout.  This helps observe logged events in realtime.
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
        os.dup2(sys.stdout.fileno(), sys.stderr.fileno())

        # execute the pipeline stage
        run_stage3(True)
