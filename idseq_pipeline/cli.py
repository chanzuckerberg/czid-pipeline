"""
idseq_pipeline

Usage:
  idseq_pipeline hello
  idseq_pipeline -h | --help
  idseq_pipeline --version

Options:
  -h --help                         Show this screen.
  --version                         Show version.

Examples:
  idseq_pipeline hello

Help:
  For help using this tool, please open an issue on the Github repository:
  https://github.com/chanzuckerberg/idseq-web
"""


from inspect import getmembers, isclass

from docopt import docopt

from . import __version__ as VERSION


def main():
    """Main CLI entrypoint."""
    import idseq_pipeline.commands
    options = docopt(__doc__, version=VERSION)

    # Here we'll try to dynamically match the command the user is trying to run
    # with a pre-defined command class we've already created.
    for (k, v) in options.items(): 
        if hasattr(idseq_pipeline.commands, k) and v:
            module = getattr(idseq_pipeline.commands, k)
            idseq_pipeline.commands = getmembers(module, isclass)
            command = [command[1] for command in idseq_pipeline.commands if command[0] != 'Base'][0]
            command = command(options)
            command.run()
