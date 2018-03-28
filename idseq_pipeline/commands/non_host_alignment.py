"""The non-host-alignment command."""
from .base import Base

class Non_host_alignment(Base):
    """Perform non-host alignment"""

    def run(self):
        from .non_host_alignment_functions import *

        unbuffer_stdout()
        upload_commit_sha(self.version)

        run_stage2(False)
