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
        """
        Class to detect atoms across a bond.
        :param in_pdb: path to pdb file
        :type in_pdb: str
        :param bond_to_descend: tuple with the atom names that define the bond
        :type bond_to_descend: tuple
        :param chain_ligand: label with the chain name (of the ligand)
        :type chain_ligand: str
        """
        pdm.PDB.__init__(self, in_pdb=in_pdb, chain=chain_ligand)
        self.ligand = self.get_ligand()
        self.bonds = self.get_bonds()
        self.names_dictionary = self.get_names_dictionary_from_ligand()
        self.bond_to_descend = bond_to_descend  # Must be a tuple

    def get_ligand(self):
        """
        Selects the ligand (by the chain) and return it.
        :return:
        """
        prody_pdb = prody.parsePDB(self.in_pdb)
        ligand = prody_pdb.select("chain {}".format(self.chain))
        return ligand

    def extract_ligand(self, path):
        """
        Writes the ligand into a file.
        :param path: path to write the ligand.
        :return:
        """
        prody.writePDB(path, self.ligand)

    def tmp_ligand(self):
        """
        Check and creates a temporal folder, and writes the ligand into it.
        :return:
        """
        if not os.path.exists("tmp"):
            create_tmp_folder()
        self.extract_ligand("tmp/ligand.pdb")

    def get_ligand_in_rdkit(self):
        """
        Reads the ligand in rdkit after writing it in a temporary folder. It is required to extract the ligand into a
        PDB file to then read it with rdkit, otherwise the reading it is not possible.
        :return: ligand as rdkit object (Mol object).
        """
        self.tmp_ligand()
        ligand = load_ligand_rdkit("tmp/ligand.pdb")
        return ligand

    def get_bonds(self):
        """
        It returns a list of atom indexes for all atoms that are doing bonds in the molecule.
        :return: list of tuples with the indexes for each bond
        """
        ligand = self.get_ligand_in_rdkit()
        bonds = get_bonds_of_molecule(ligand)
        return bonds

    def check_if_bond_exist(self, bond):
        """
        It checks if the selected bond exist in the molecule.
        :param bond: bond (tuple with atom's indexes)
        :return: checks if exist, otherwise raises an error.
        """
        if bond not in self.bonds and tuple(reversed(bond)) not in self.bonds:
            raise ValueError("Bond {} does not exist!".format(bond))

    def get_descendent_tree_from_bond(self, bond_to_descend_from):
        """
        Transforms the molecule into a graph, where edges are bonds and nodes are atoms. It visit all the nodes (atoms)
        that goes after the bond_to_descend_from (initial edge), and returns them.
        :param bond_to_descend_from: bond (tuple with atom's indexes)
        :return: list of nodes or atoms that goes after the bond_to_descend_from
        """
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
        """
        Given an atom (node), it computes all possible atoms (nodes) that must be visited (the ones that have not been
        visited yet).
        :param atom: index of an atom
        :param atoms_visited: list of indexes with the atoms that have been visited
        :param atoms_to_visit: list of indexes with the atoms to visit
        :return: updates atoms_visited and atoms_to visit
        """
        for bond in self.bonds:
            if atom in bond and atom not in atoms_visited:  # If the atoms has not been visited...
                for bond_atom in bond:
                    if bond_atom != atom:
                        atoms_to_visit.append(bond_atom)  # Put it into the list of atoms to visit
        atoms_visited.add(atom)  # If its here is because it has been visited, so we update the list
        return atoms_visited, atoms_to_visit

    def get_names_dictionary_from_ligand(self):
        if not os.path.exists("tmp/ligand.pdb"):
            self.tmp_ligand()
        pdb_ligand = pdm.PDB("tmp/ligand.pdb", self.chain)
        dictionary_names = pdb_ligand.get_names_dictionary()
        return dictionary_names

    def bond_names_to_indexes(self):
        """
        Transform the pdb atom names of the bond_to_descend to atom indexes
        :return: tuple of indexes
        """
        dictionary_of_names = self.get_names_dictionary_from_ligand()
        indexes = [None, None]
        for index, name in dictionary_of_names.items():
            if name == self.bond_to_descend[0]:
                indexes[0] = int(index)
            if name == self.bond_to_descend[1]:
                indexes[1] = int(index)
        return tuple(indexes)


def get_bonds_of_molecule(molecule):
    """
    It returns the atom's indexes for all bonds of a molecule.
    :param molecule: Mol object with the ligand (Rdkit)
    :return: list of tuples with the indexes for each bond
    """
    bonding_idx = []
    for bond in molecule.GetBonds():
        bonding_idx.append((bond.GetEndAtomIdx() + 1, bond.GetBeginAtomIdx() + 1))
    return bonding_idx


def load_ligand_rdkit(path_to_ligand):
    """
    Loads the ligand as Mol object. It keeps the hydrogen atoms and don't sanitize because we want to keep the molecule
    as it is in the PDB.
    :param path_to_ligand: path to the ligand pdb
    :type path_to_ligand: str
    :return: Mol object with the ligand
    """
    mol = Chem.MolFromPDBFile(path_to_ligand, removeHs=False, sanitize=False)
    return mol


def create_tmp_folder():
    os.mkdir("tmp")


def remove_tmp_folder():
    shutil.rmtree("tmp")


def main(pdb_complex, bond_to_descend, chain_ligand="L"):
    """
    Returns the pdb atom names of the atoms that goes after the selected bond (bond_to_descend).
    :param pdb_complex: path to the pdb complex
    :type pdb_complex: str
    :param bond_to_descend: tuple with the pdb atom names that define the selected bond
    :type bond_to_descend: tuple
    :param chain_ligand: label of the chain with the ligand
    :type chain_ligand: str
    :return: list of pdb atom names after the selected bond
    """
    atoms_in_tree = []
    to_detect = Detector(in_pdb=pdb_complex, bond_to_descend=bond_to_descend, chain_ligand=chain_ligand)
    bond_to_descend_indexes = to_detect.bond_names_to_indexes()
    bonds_detected_indexes = to_detect.get_descendent_tree_from_bond(bond_to_descend_indexes)
    for index in bonds_detected_indexes:
        name = to_detect.names_dictionary[index]
        atoms_in_tree.append(name)
    remove_tmp_folder()
    return atoms_in_tree
