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


Developers
----------


When merging a commit to master, you need to increase the version number in `idseq_pipeline/version/__init__.py`:
  - if results are expected to change, increase the 2nd number
  - if results are not expected to change, increase the 3rd number.


Release notes
-------------

- 1.6.0
  - Update to GMAP-GSNAP 2018-03-25.

- 1.5.0
    - Identify hits that match multiple species within the same genus as
      "non species specific" hits to the genus.

- 1.4.0
    - Version result folder.

- 1.3.0
    - Fix bug causing alignment to run before host subtraction in samples
      with unpaired reads.
    - Include ERCC gene counts from STAR.

- 1.2.0
    - Synchronize pair order after STAR to improve sensitivity in 10% of
      samples with paired-end reads.
