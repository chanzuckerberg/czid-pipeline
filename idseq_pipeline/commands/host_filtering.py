"""The host-filtering command."""
from .base import Base
import sys

class Host_filtering(Base):
    """Perform host filtering"""

    def run(self):
        from .host_filtering_functions import *

        unbuffer_stdout()
        upload_commit_sha(self.version, SAMPLE_S3_OUTPUT_PATH)
        upload_pipeline_version_file()

        print("Python version")
        print(sys.version)
        print("Version info.")
        print(sys.version_info)

        run_stage1(True)
