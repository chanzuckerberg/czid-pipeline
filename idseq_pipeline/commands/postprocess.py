"""The post-processing command."""
from .base import Base


class Postprocess(Base):
    """Perform post-processing"""

    def run(self):
        from .postprocess_functions import *

        set_up_stdout()
        set_up_commit_sha(self.version)

        # TODO(yf): didn't do any pipeline version check. Revisit later

        try:
          run_stage3(lazy_run=False)
          mark_job_complete(SAMPLE_S3_OUTPUT_PATH, JOB_SUCCEEDED)
        except Exception as e:
          mark_job_complete(SAMPLE_S3_OUTPUT_PATH, JOB_FAILED)
          raise RuntimeError("Job failed: ", e)
