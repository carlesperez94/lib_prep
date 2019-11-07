class PDB:
    def __init__(self, in_pdb, chain="L", resname="LIG", resnum="   1"):
        self.in_pdb = in_pdb
        self.content = self.read_content_as_str()
        self.lines = self.read_content_as_lines()
        self.atom_section = self.read_atoms_section_as_str()
        self.conect_section = self.read_conect()
        self.chain = chain
        self.resname = resname
        self.resnum = resnum

    def read_content_as_str(self):
        with open(self.in_pdb) as infile:
            content = infile.read()
        return content

    def read_content_as_lines(self):
        with open(self.in_pdb) as infile:
            content = infile.readlines()
        return content

    def read_atoms_section(self):
        atoms_sect = []
        for line in self.lines:
            if "ATOM" in line[0:4] or "HETATM" in line[0:6]:
                atoms_sect.append(line)
        return atoms_sect

    def read_conect(self):
        conect_sect = []
        for line in self.lines:
            if "CONECT" in line[0:6]:
                conect_sect.append(line)
        return conect_sect

    def read_atoms_section_as_str(self):
        return "\n".join(self.read_atoms_section())

    def modify_pdb(self):
        pdb_atom_lines = self.read_atoms_section()
        new_pdb = []
        conect_lines = self.conect_section
        for line in pdb_atom_lines:
            line2 = set_resname_to_line(line, self.resname)
            line3 = set_resnum_to_line(line2, self.resnum)
            line4 = set_chain_to_line(line3, self.chain)
            new_pdb.append(line4)
        self.content = "".join(new_pdb)+"".join(conect_lines)
        return self.content

    def overwrite_pdb(self):
        with open(self.in_pdb, "w") as out_pdb:
            out_pdb.write(self.content)
            print("{} has been modified.".format(self.in_pdb))

    def find_pdb_atom_name_of_idx(self, index):
        for line in self.read_atoms_section():
            if int(get_atom_index_from_line(line)) == int(index):
                return get_atom_pdb_name_from_line(line)

    def get_names_dictionary(self):
        names_dict = {}
        for line in self.read_atoms_section():
            index = int(get_atom_index_from_line(line).strip())
            name = get_atom_pdb_name_from_line(line).strip()
            names_dict[index] = name
        return names_dict

    def get_atoms_of_chain(self):
        atoms_of_chain = []
        for line in self.read_atoms_section():
            chain = get_chain_from_line(line).strip()
            if chain == self.chain:
                atoms_of_chain.append(line)
        return atoms_of_chain


def get_resname_from_line(line):
    return line[17:20]


def get_chain_from_line(line):
    return line[21:22]


def get_resnum_from_line(line):
    return line[22:26]


def get_atom_index_from_line(line):
    return line[6:11]


def get_atom_pdb_name_from_line(line):
    return line[12:16]


def set_resname_to_line(line, value):
    line = list(line)
    line[17:20] = value
    line = "".join(line)
    return line


def set_chain_to_line(line, value):
    line = list(line)
    line[21:22] = value
    line = "".join(line)
    return line


def set_resnum_to_line(line, value):
    line = list(line)
    line[22:26] = "{:>4}".format(value)
    line = "".join(line)
    return line


def set_index_to_line(line, value):
    line = list(line)
    line[7:11] = "{:>4}".format(value)
    line = "".join(line)
    return line







