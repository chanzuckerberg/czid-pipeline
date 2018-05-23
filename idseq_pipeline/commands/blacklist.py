"""The blacklist command."""
from .base import Base


class Blacklist(Base):
    """Generate blacklist"""

    def run(self):
        from .blacklist_functions import *

        set_up_stdout()

        # make blacklist
        make_blacklist()
