"""The non-host-alignment command."""
from .base import Base

class Non_host_alignment(Base):
    """Perform non-host alignment"""

    def run(self):
        from .non_host_alignment_functions import *

        unbuffer_stdout()
        lazy_run = not big_version_change_from_last_run(self.version, get_alignment_version_s3_path())
        upload_commit_sha(self.version)

        run_stage2(lazy_run)
