"""The host-filtering command."""
from .base import Base

class Host_filtering(Base):
    """Perform host filtering"""

    def run(self):
        from .host_filtering_functions import *

        unbuffer_stdout()
        lazy_run = not big_version_change_from_last_run(self.version, get_host_filtering_version_s3_path())
        upload_commit_sha(self.version)

        run_stage1(lazy_run)
