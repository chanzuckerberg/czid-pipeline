"""The blacklist command."""
from .base import Base

class Blacklist(Base):
    """Generate blacklist"""

    def run(self):
        from .blacklist_functions import *

        # Unbuffer stdout and redirect stderr into stdout.  This helps observe logged events in realtime.
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
        os.dup2(sys.stdout.fileno(), sys.stderr.fileno())

        # make blacklist
        make_blacklist()
