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

- 1.6.1
    - Perform de-novo assembly using SPAdes for species with >= 100 reads.

- 1.6.0
    - Fix fasta downloads broken by release 1.5.0, making sure only
      hits at the correct level are output in the deduped m8.
    - Fix fasta download for samples with unpaired reads by eliminating
      merged fasta for those samples.
    - Extend the partial fix in release 1.5.1 to repair more of the
      broken reports.  Full fix requires rerun with updated webapp.
    - Correctly aggregate counts for species with unclassified genera,
      such as e.g. genus-less species 1768803 from family 80864.
    - Fix total count in samples with unpaired reads (no longer doubled).
    - Fix crash when zero reads remain after host filtering.
    - Fix bug in enforcing command timeouts that could lead to hangs.
    - Fix performance regression in stage 2 (non-host alignment)
      introduced with 1.5.0.
    - Deduplicate and simplify much of stage 2, and improve performance
      by parallelizing uploads and downloads.

- 1.5.1
    - Fix bug introduced in 1.5.0 breaking samples with non-species-specific
      deuterostome hits.

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
