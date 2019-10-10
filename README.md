# LibPrep
Python package to analyse and prepare libraries of chemical compounds for molecular simulations.

## Dependences
- Python 3.7
- Rdkit 2019.03.4
- Schrodinger 2019 (Optional)

# Programs available
## convert_sdfs2pdb.py
Program to transform 3D SDF files with large amount of molecular compounds to several PDBs, renaming also chains, resnames, and resnums to the desired ones.

Example of usage:

python convert_sdfs2pdb.py library.sdf --out_folder library_in_pdb --chain L --resname LIG --resnum 1

For more information:

python convert_sdfs2pdb.py -h

## FragmentTools/prepare_to_frag.py
Program to prepare FragPELE instruction's from PDB libraries.

Exaple of usage:

python FragmentTools/prepare_to_frag.py /path/to/library /path/to/complex.pdb C1 -m first-occurrence -o /out/path/instructions.conf

For more information:

python FragmentTools/prepare_to_frag.py -h
