==============
Documentation
==============

.. toctree::
   :maxdepth: 2


Convert2PDB
-----------

General usage::

   convert_sdfs2pdb.py [-h] [-o OUT_FOLDER] [-ch CHAIN] [-res RESNAME] [-rn RESNUM] [-pn NAMING_PROPERTY] sdf_input

Required parameters:
++++++++++++++++++++
Take in mind that these parameters are positionals, so the order matters:

- **sdf_input**: Path to a library of fragments or ligands in SDF file. The ligands/fragments must be previously protonated and optimized .

Optional arguments:
+++++++++++++++++++
-  **-o**, **---out_folder**: Path of the output folder. PDB output files will be stored in this path. Default: "out". Type: str.
-  **-ch**, **---chain**: Label for the chain of the ligand in new PDBs. Default: "L". Type: str.
-  **-res**, **---resname**: Label for the residue name of the ligand in new PDBs. Default: "LIG". Type: str.
-  **-rn**, **---resnum**: Label for the residue number of the ligand in new PDBs. Default: 1. Type: int.
-  **-pn**, **---naming_property**: Property of the SDF file that will we used to name the output PDB files. By default: "Molecule Name". If the property is not found in the SDF file the molecules are renamed automatically as "MOL{:09d}", where the content in brackets will be replaced by a 9 digits counter.

Example:
++++++++

Normally you do not need to change any optional parameter to run Convert2PDB, so for common use try the following command::

   convert_sdfs2pdb.py library.sdf

If, for some reason, you want that to name the chain of the ligands as "Z", the residue name as "WAT", and the residue number as 100::

   convert_sdfs2pdb.py library.sdf --chain Z --resname WAT --resnum 100

Prepare To Frag
---------------

General usage::

   prepare_to_frag.py [-h] [-m {first-occurrence}] [-o OUT] library_path pdb_complex heavy_atom_pdb_complex

Required parameters:
++++++++++++++++++++
Take in mind that these parameters are positionals, so the order matters:

- **library_path**: Path to the folder that contains the library of fragments or ligands in PDF file. This ideally would be the output folder of Convert2PDB, but if you have another library of PDB files already prepared it can be used anyway. REMINDER: the program will check if all files inside the folder are PDBs, otherwise an exception will rise!

- **pdb_complex**: Path to a ligand-protein complex. The ligand should be the core or scaffold to grow fragments onto it with FragPELE.

- **heavy_atom_pdb_complex**: Atom name of the heavy atom of the previous core/scaffold to grow fragments onto it. All links between core and fragment will be generated through this atom.

Optional arguments:
+++++++++++++++++++

-  **-o**, **---out**: Path of the output file. By default the output file will be stored into the library_path and named "serie_file_{}.conf", filling the brackets with the pdb_complex name. Type: str.
-  **-m**, **---mode**: Choose between different scanning modes.
   - first-occurrence: the first heavy atom found with at least one hydrogen atom bonded will be selected.

Currently, only one single mode is available but in the future we will implement more.

Example:
++++++++

Easy example of usage::

   FragmentTools/prepare_to_frag.py /out/pdbs /path/to/complex.pdb C1
