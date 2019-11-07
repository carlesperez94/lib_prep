#!/usr/bin/env python
"""
Program that given a protein-ligand complex (PDB) and a bond of the ligand (tuple with PDB ATOM NAMES), it returns the
PDB ATOM NAMES of the atoms that are across the bond.
"""

import os
import shutil
import prody
from rdkit import Chem
import lib_prep.pdb_modifier as pdm

__author__ = "Carles Perez Lopez"
__version__ = "1.0.0"
__maintainer__ = "Carles Perez Lopez"
__email__ = "carlesperez94@gmail.com"


class Detector(pdm.PDB):
    def __init__(self, in_pdb, bond_to_descend, chain_ligand="L"):
        pdm.PDB.__init__(self, in_pdb=in_pdb, chain=chain_ligand)
        self.ligand = self.get_ligand()
        self.bonds = self.get_bonds()
        self.names_dictionary = self.get_names_dictionary_from_ligand()
        self.bond_to_descend = bond_to_descend  # Must be a tuple

    def get_ligand(self):
        prody_pdb = prody.parsePDB(self.in_pdb)
        ligand = prody_pdb.select("chain {}".format(self.chain))
        return ligand

    def extract_ligand(self, path):
        prody.writePDB(path, self.ligand)

    def tmp_ligand(self):
        if not os.path.exists("tmp"):
            create_tmp_folder()
        self.extract_ligand("tmp/ligand.pdb")

    def get_ligand_in_rdkit(self):
        self.tmp_ligand()
        ligand = load_ligand_rdkit("tmp/ligand.pdb")
        return ligand

    def get_bonds(self):
        ligand = self.get_ligand_in_rdkit()
        bonds = get_bonds_of_molecule(ligand)
        return bonds

    def check_if_bond_exist(self, bond):
        if bond not in self.bonds and tuple(reversed(bond)) not in self.bonds:
            raise ValueError("Bond {} does not exist!".format(bond))

    def get_descendent_tree_from_bond(self, bond_to_descend_from):  # O - C
        self.check_if_bond_exist(bond_to_descend_from)
        atoms_visited = set()
        atoms_visited.add(bond_to_descend_from[0])  # This node has to be added to avoid double direction in the tree
        atoms_to_visit = [bond_to_descend_from[1]]  # This is the first node that should be visited
        while atoms_to_visit:
            atom = atoms_to_visit.pop()
            atoms_visited, atoms_to_visit = self.get_branches_of_tree_from_bond(atom, atoms_visited, atoms_to_visit)
        atoms_visited.remove(bond_to_descend_from[0])
        return list(atoms_visited)

    def get_branches_of_tree_from_bond(self, atom, atoms_visited, atoms_to_visit):
        for bond in self.bonds:
            if atom in bond and atom not in atoms_visited:
                for bond_atom in bond:
                    if bond_atom != atom:
                        atoms_to_visit.append(bond_atom)
        atoms_visited.add(atom)
        return atoms_visited, atoms_to_visit

    def get_names_dictionary_from_ligand(self):
        if not os.path.exists("tmp/ligand.pdb"):
            self.tmp_ligand()
        pdb_ligand = pdm.PDB("tmp/ligand.pdb", self.chain)
        dictionary_names = pdb_ligand.get_names_dictionary()
        return dictionary_names

    def bond_names_to_indexes(self):
        dictionary_of_names = self.get_names_dictionary_from_ligand()
        indexes = [None, None]
        for index, name in dictionary_of_names.items():
            if name == self.bond_to_descend[0]:
                indexes[0] = int(index)
            if name == self.bond_to_descend[1]:
                indexes[1] = int(index)
        return tuple(indexes)


def get_bonds_of_molecule(molecule):
    bonding_idx = []
    for bond in molecule.GetBonds():
        bonding_idx.append((bond.GetEndAtomIdx() + 1, bond.GetBeginAtomIdx() + 1))
    return bonding_idx


def load_ligand_rdkit(path_to_ligand):
    mol = Chem.MolFromPDBFile(path_to_ligand, removeHs=False, sanitize=False)
    return mol


def create_tmp_folder():
    os.mkdir("tmp")


def remove_tmp_folder():
    shutil.rmtree("tmp")


def main(pdb_complex, bond_to_descend, chain_ligand="L"):
    atoms_in_tree = []
    to_detect = Detector(in_pdb=pdb_complex, bond_to_descend=bond_to_descend, chain_ligand=chain_ligand)
    bond_to_descend_indexes = to_detect.bond_names_to_indexes()
    bonds_detected_indexes = to_detect.get_descendent_tree_from_bond(bond_to_descend_indexes)
    for index in bonds_detected_indexes:
        name = to_detect.names_dictionary[index]
        atoms_in_tree.append(name)
    remove_tmp_folder()
    return atoms_in_tree
