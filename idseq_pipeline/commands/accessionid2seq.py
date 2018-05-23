"""The  accessionid2seq command."""
from .base import Base


class Accessionid2seq(Base):
    """Accessioid2seq"""

    def run(self):
        from .accessionid2seq_functions import accessionid2seq_main
        from .common import *

        set_up_stdout()

        accessionid2seq_main(arguments=self.options)
