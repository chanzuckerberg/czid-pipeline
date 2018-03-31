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


1.3.0 -- Correctly subtract host for non-paired reads.  Previously for
         non-paired reads, alignment was being performed on a random
         sample of the unfiltered input.

         Include ERCC gene counts from STAR.

1.2.0 -- Synchronize misordered pairs after STAR.
