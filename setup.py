"""Packaging settings."""


from codecs import open
from os.path import abspath, dirname, join
from subprocess import call

from setuptools import Command, find_packages, setup

from idseq_pipeline import __version__


this_dir = abspath(dirname(__file__))
with open(join(this_dir, 'README.rst'), encoding='utf-8') as file:
    long_description = file.read()


class RunTests(Command):
    """Run all tests."""
    description = 'run tests'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        """Run all tests!"""
        errno = call(['py.test', '--cov=idseq_pipeline', '--cov-report=term-missing'])
        raise SystemExit(errno)


setup(
    name = 'idseq_pipeline',
    version = __version__,
    description = 'IDseq pipeline',
    long_description = long_description,
    url = 'https://github.com/chanzuckerberg/idseq-pipeline',
    author = 'idseq-pipeline contributors',
    author_email = 'cdebourcy@chanzuckerberg.com',
    license = open("LICENSE").readline().strip(),
    classifiers = [
        'Intended Audience :: Developers',
        'Topic :: Utilities',
        'License :: Public Domain',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    keywords = 'cli',
    packages = find_packages(exclude=['docs', 'tests*']),
    install_requires = [line.rstrip() for line in open(os.path.join(os.path.dirname(__file__), "requirements.txt"))],
    extras_require = {
        'test': ['coverage', 'pytest', 'pytest-cov'],
    },
    entry_points = {
        'console_scripts': [
            'idseq_pipeline=idseq_pipeline.cli:main',
        ],
    },
    cmdclass = {'test': RunTests},
)
