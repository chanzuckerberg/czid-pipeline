"""
idseq_pipeline
 
Usage:
  idseq_pipeline host_filtering
  idseq_pipeline non_host_alignment
  idseq_pipeline postprocess
  idseq_pipeline indexing
  idseq_pipeline blacklist
  idseq_pipeline -h | --help
  idseq_pipeline --version
 
Options:
  -h --help                         Show this screen.
  --version                         Show version.
 
Examples:
  idseq_pipeline host_filtering
 
Help:
  For help using this tool, please open an issue on the Github repository:
  https://github.com/chanzuckerberg/idseq-pipeline
"""
 
 
from inspect import getmembers, isclass
 
from docopt import docopt
 
from . import __version__ as VERSION
 
 
def main():
    """Main CLI entrypoint."""
    import commands
    options = docopt(__doc__, version=VERSION)
 
    # Here we'll try to dynamically match the command the user is trying to run
    # with a pre-defined command class we've already created.
    for k, v in options.iteritems():
        if hasattr(commands, k):
            module = getattr(commands, k)
            commands = getmembers(module, isclass)
            command = [command[1] for command in commands if command[0] != 'Base'][0]
            command = command(options)
            command.run()
