===============
Getting Started
===============

.. toctree::
   :maxdepth: 2


The general workflow of the preparation of a fragment's library follows the next schema:

.. image:: img/schema.png
   :scale: 80 %
   :align: center

Input Files
-----------

Therefore, the required input file is:

* SDF file: Must contain fragments or ligands in 3D, optimized, and protonated. We recommend to use previously LigPrep onto the library to do the preparation.


Pipeline
--------
The whole pipeline consist of two main steps: prepare and scan. The first one would be performed by the script convert2pdb.py,
and the second by FragmentTools/prepare2frag.py

**Convert2PDB**
+++++++++++++++

Program to transform 3D SDF files with large amount of molecular compounds to several PDBs, renaming also chains,
resnames, and resnums to the desired ones.

i.e::

    python convert_sdfs2pdb.py library.sdf --out_folder library_in_pdb --chain L --resname LIG --resnum 1


**Prepare to Frag**
+++++++++++++++++++

Program to prepare FragPELE instruction's from PDB libraries.

i.e::

    python FragmentTools/prepare_to_frag.py /path/to/library /path/to/complex.pdb C1 -m first-occurrence -o /out/path/instructions.conf

