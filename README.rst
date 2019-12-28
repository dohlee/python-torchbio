========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |appveyor| |requires|
        | |codecov|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/python-torchbio/badge/?style=flat
    :target: https://readthedocs.org/projects/python-torchbio
    :alt: Documentation Status

.. |travis| image:: https://api.travis-ci.org/dohlee/python-torchbio.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/dohlee/python-torchbio

.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/github/dohlee/python-torchbio?branch=master&svg=true
    :alt: AppVeyor Build Status
    :target: https://ci.appveyor.com/project/dohlee/python-torchbio

.. |requires| image:: https://requires.io/github/dohlee/python-torchbio/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/dohlee/python-torchbio/requirements/?branch=master

.. |codecov| image:: https://codecov.io/github/dohlee/python-torchbio/coverage.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/dohlee/python-torchbio

.. |version| image:: https://img.shields.io/pypi/v/torchbio.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/torchbio

.. |wheel| image:: https://img.shields.io/pypi/wheel/torchbio.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/torchbio

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/torchbio.svg
    :alt: Supported versions
    :target: https://pypi.org/project/torchbio

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/torchbio.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/torchbio

.. |commits-since| image:: https://img.shields.io/github/commits-since/dohlee/python-torchbio/v0.0.0.svg
    :alt: Commits since latest release
    :target: https://github.com/dohlee/python-torchbio/compare/v0.0.0...master



.. end-badges

Pytorch utilities for bioinformatics.

* Free software: MIT license

Installation
============

::

    pip install torchbio

You can also install the in-development version with::

    pip install https://github.com/dohlee/python-torchbio/archive/master.zip


Documentation
=============


https://python-torchbio.readthedocs.io/


Development
===========

To run the all tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
