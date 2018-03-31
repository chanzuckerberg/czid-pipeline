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


Release notes
-------------

- 1.3.0   
    - Fix bug causing alignment to run before host subtraction in samples
      with unpaired reads.
    - Include ERCC gene counts from STAR.

- 1.2.0
    - Work around STAR concurrency bug causing misordered pairs in 10%
      of samples with paired end reads, causing reduced sensitivity
      in gsnap etc.
