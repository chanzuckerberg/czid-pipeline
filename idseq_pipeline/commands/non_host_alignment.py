#!/usr/bin/env python3

"""The non-host-alignment command."""
from .base import Base


class Non_host_alignment(Base):
    """Perform non-host alignment"""

    def run(self):
        from .non_host_alignment_functions import *
        unbuffer_stdout()
        upload_commit_sha(self.version)
        stage1_version = os.environ.get('HOST_FILTER_PIPELINE_VERSION')
        if stage1_version and major_version(
                self.version) != major_version(stage1_version):
            raise RuntimeError(
                "Stage 1 and 2 run on different version: %s vs %s" %
                (stage1_version, self.version))

        print("Python version")
        print(sys.version)
        print("Version info.")
        print(sys.version_info)

        run_stage2(True)
