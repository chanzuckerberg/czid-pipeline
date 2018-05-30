"""The host-filtering command."""
from .base import Base


class Host_filtering(Base):
    """Perform host filtering"""

    def run(self):
        from .host_filtering_functions import *

        set_up_stdout()
        set_up_commit_sha(self.version, SAMPLE_S3_OUTPUT_PATH)
        upload_pipeline_version_file()

        try:
          run_stage1(True)
          mark_job_complete(SAMPLE_S3_OUTPUT_PATH, JOB_SUCCEEDED)
        except:
          mark_job_complete(SAMPLE_S3_OUTPUT_PATH, JOB_FAILED)
          raise RuntimeError("Job failed")

