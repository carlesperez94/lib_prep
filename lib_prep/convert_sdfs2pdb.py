#!/usr/bin/env python
"""
Program to transform 3D SDF files with large amount
of molecular compounds to several PDBs.
"""
import os
import argparse
from lib_prep.LibraryManager import lib_manager, library


__author__ = "Carles Perez Lopez"
__version__ = "1.0.0"
__maintainer__ = "Carles Perez Lopez"
__email__ = "carlesperez94@gmail.com"


def parse_arguments():
    """
        Parse user arguments
        Output: list with all the user arguments
    """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""Description: Program to transform 3D SDF files with large amount
    of molecular compounds to several PDBs.""")
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument("sdf_input", type=str, help="SDF input file. Must contain 3D data.")
    parser.add_argument("-o", "--out_folder", type=str, default="out",
                        help="Path to the output folder.")
    parser.add_argument("-ch", "--chain", type=str, default="L",
                        help="Chain name for all pdbs in the library. Max char length: 1.")
    parser.add_argument("-res", "--resname", type=str, default="LIG",
                        help="Resname for all pdbs in the library. Max char length: 4")
    parser.add_argument("-rn", "--resnum", type=str, default="   1",
                        help="Resnum for all pdbs in the library. Max char length: 4")
    parser.add_argument("-pn", "--naming_property", type=str, default="Molecule Name",
                        help="Name of the property of the SDF file that will be used to rename the output PDB files.")

    args = parser.parse_args()

    return args.sdf_input, args.out_folder, args.chain, args.resname, args.resnum, args.naming_property


def main(sdf_input, output_folder, chain, resname, resnum, property_to_name):
    lib = library.LibrarySDF(sdf_input)
    lib.export_sdf_to_pdb_rdkit(out_path=output_folder, molname_property=property_to_name)
    pdbs_library = os.path.join(output_folder, "pdbs")
    pdb_lib = lib_manager.LibraryChecker(pdbs_library)
    pdb_lib.prepare_library(chain, resname, resnum)


if __name__ == '__main__':
    sdf_in, output_fold, chain, resname, resnum, property_to_name= parse_arguments()
    main(sdf_input=sdf_in, output_folder=output_fold, chain=chain, resname=resname, resnum=resnum,
         property_to_name=property_to_name)

