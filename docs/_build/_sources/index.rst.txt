.. LibPrep documentation master file, created by
   sphinx-quickstart on Tue Nov 12 16:22:00 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

LibPrep: Library Preparer
==========================

LibPrep is a Python Package with several modules to process PDB files (complexes) and SDF files (drug libraries).

Fragment's libraries usually are huge files that contains thousands of chemical compounds, with variable chemical
properties. LibPrep aims to facilitate the processing and pre-analysis of all these molecules and filter them to
further uses, such as running more computationally expensive programs (`PELE <https://pele.bsc.es/pele.wt>`_).

It has been designed to automatically prepare the input files to run `FragPELE <https://github.com/carlesperez94/frag_pele/tree/master>`_
towards a fragment library in SDF format. However, it can be used with other proposes.

The general workflow of the preparation of a fragment's library follows the next schema:

.. image:: img/schema.png
   :scale: 80 %
   :align: center

LibPrep is focused in the "Prepare" and "Scan" parts of the process. In previous steps, you will require LigPrep or
another similar program to protonate and optimize the molecules to the correct pH and 3D geometry.

Requirements:

- Python 3.6 or higher
- Pandas 0.18.0 or higher
- Rdkit 2019.03.4 or higher

Installation & Starting
=======================
.. toctree::
   installation
   first_steps

Documentation
=============
.. toctree::
   documentation



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
