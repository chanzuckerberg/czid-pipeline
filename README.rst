idseq_pipeline
=========

*The data processing pipeline for IDseq.*


Purpose
-------

This is a CLI that allows you to execute the different data processing stages required in IDseq.


Usage
-----

To install after cloning:

    $ pip install -e .

When merging a commit to master, you need to increase the version number in `idseq_pipeline/version/__init__.py`:
  - if results are expected to change, increase the 2nd number
  - if results are not expected to change, increase the 3rd number.
