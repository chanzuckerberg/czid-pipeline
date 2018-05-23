"""The non-host-alignment command."""
from .base import Base


class Non_host_alignment(Base):
    """Perform non-host alignment"""

    def run(self):
        from .non_host_alignment_functions import *
        set_up_stdout()
        set_up_commit_sha(self.version)
        stage1_version = os.environ.get('HOST_FILTER_PIPELINE_VERSION')
        if stage1_version and major_version(
                self.version) != major_version(stage1_version):
            raise RuntimeError(
                "Stage 1 and 2 run on different version: %s vs %s" %
                (stage1_version, self.version))

        run_stage2(True)
