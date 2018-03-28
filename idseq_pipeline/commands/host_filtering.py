"""The host-filtering command."""
from .base import Base

class Host_filtering(Base):
    """Perform host filtering"""

    def run(self):
        from .host_filtering_functions import *

        unbuffer_stdout()
        upload_commit_sha(self.version)

        run_stage1(False)
