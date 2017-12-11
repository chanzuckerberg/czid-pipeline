import os
from setuptools import setup

from idseq_pipeline import __version__

install_requires = [line.rstrip() for line in open(os.path.join(os.path.dirname(__file__), "requirements.txt"))]

setup(
    name="idseq_pipeline",
    version=__version__,
    url='https://github.com/chanzuckerberg/idseq-pipeline',
    license=open("LICENSE").readline().strip(),
    author='idseq-pipeline contributors',
    author_email='cdebourcy@chanzuckerberg.com',
    description='IDseq pipeline',
    long_description=open('README.md').read(),
    install_requires=install_requires,
    extras_require={},
    packages=['idseq_pipeline'],
    zip_safe=False,
    entry_points={'console_scripts': ['idseq_pipeline=idseq_pipeline.cli:main']}
)
