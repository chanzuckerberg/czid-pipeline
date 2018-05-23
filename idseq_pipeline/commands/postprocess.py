"""The post-processing command."""
from .base import Base


class Postprocess(Base):
    """Perform post-processing"""

    def run(self):
        from .postprocess_functions import *

        set_up_stdout()
        set_up_commit_sha(self.version)

        # TODO(yf): didn't do any pipeline version check. Revisit later

        run_stage3(lazy_run=False)
