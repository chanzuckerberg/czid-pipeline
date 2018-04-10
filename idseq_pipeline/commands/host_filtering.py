"""The host-filtering command."""
from .base import Base

class Host_filtering(Base):
    """Perform host filtering"""

    def run(self):
        from .host_filtering_functions import *

        unbuffer_stdout()
        upload_commit_sha(self.version, SAMPLE_S3_OUTPUT_PATH)
        upload_pipeline_version_file()

        run_stage1(False)
