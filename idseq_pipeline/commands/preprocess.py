"""Pre-processing commands."""
from .base import Base

class Preprocess(Base):
    """Perform pre-processing functions"""

    def run(self):
        from .preprocess_functions import *

        unbuffer_stdout()
