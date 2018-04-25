"""The assembly command."""
from .base import Base

class Assembly(Base):
    """Perform assembly"""

    def run(self):
        from .assembly_functions import *

        unbuffer_stdout()
        upload_commit_sha(self.version)

        run_stage4()
