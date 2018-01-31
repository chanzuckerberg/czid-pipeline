"""The post-processing command."""
from .base import Base

class Postprocess(Base):
    """Perform post-processing"""

    def run(self):
        from .postprocess_functions import *

        unbuffer_stdout()
        upload_commit_sha(self.version)

        run_stage3(True)
