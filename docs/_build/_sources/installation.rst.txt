.. LibPrep documentation master file, created by
   sphinx-quickstart on Tue Nov 12 16:22:00 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
==========================

.. toctree::
   :maxdepth: 2

Conda (Recommended)
-------------------

Conda installation commands::

    conda create -n libprep python=3.7
    conda activate libprep
    conda install -c NostrumBioDiscovery

Once you have the conda environment created, you just need to activate it to use lib_prep::
    
    conda activate libprep

Pypi
----

TODO

Source code
-----------

TODO

Test
----

Test LibPrep through pytest module::
    cd tests/
    python -m pytest test_all.py
