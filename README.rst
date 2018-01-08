idseq_pipeline
=========

*The data processing pipeline for IDseq.*


Purpose
-------

This is a CLI that allows you to execute the different data processing stages required in IDseq.


Usage
-----

To install after cloning:

    $ pip install -e .[test]

To run tests:

    $ python setup.py test

To cut a new release and publish it to the Python Package Index:

    $ python setup.py sdist bdist_wheel
    $ twine upload dist/*
