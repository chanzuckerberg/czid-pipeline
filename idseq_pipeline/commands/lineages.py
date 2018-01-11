"""The lineage command."""
from .base import Base

class Lineages(Base):
    """Make lineage file"""

    def run(self):
        from .common import *

        unbuffer_stdout()

        execute_command("git clone https://github.com/chanzuckerberg/ncbitax2lin.git")
        execute_command("cd ncbitax2lin; make")
        execute_command("aws s3 cp taxid-lineages.db %s/" % OUTPUT_PATH
