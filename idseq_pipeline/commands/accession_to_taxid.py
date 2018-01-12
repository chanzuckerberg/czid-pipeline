"""The accession-to-taxid mapping command."""
from .base import Base

class Accession_to_taxid(Base):
    """Make accession-to-taxid DB"""

    def run(self):
        from .accession_to_taxid_functions import *

        unbuffer_stdout()
        make_accession_db()
