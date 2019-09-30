#!/usr/bin/env python
"""
Program to transform 3D SDF files with large amount
of molecular compounds to several PDBs.
"""

import library
import argparse

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
    args = parser.parse_args()

    return args.sdf_input, args.out_folder


def main(sdf_input, output_folder):
    lib = library.Library(sdf_input)
    lib.export_sdf_to_pdbs(out_path=output_folder)


if __name__ == '__main__':
    sdf_in, output_fold = parse_arguments()
    main(sdf_input=sdf_in, output_folder=output_fold)

