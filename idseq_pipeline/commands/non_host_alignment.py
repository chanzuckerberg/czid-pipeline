"""The non-host-alignment command."""
from .base import Base

class Non_host_alignment(Base):
    """Perform non-host alignment"""

    def run(self):
        from .non_host_alignment_functions import *
        unbuffer_stdout()
        upload_commit_sha(self.version)
        stage1_major_version = os.environ.get('HOST_FILTER_PIPELINE_VERSION')
        if stage1_major_version and major_version(self.version) != stage1_major_version:
            raise Exception("Stage 1 and 2 run on different version: %s vs %s" % (stage1_major_version, self.version))

        run_stage2(True)
