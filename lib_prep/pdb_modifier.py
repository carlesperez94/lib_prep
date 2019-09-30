class PDBModifier:
    def __init__(self, in_pdb):
        self.in_pdb = in_pdb
        self.content = self.read_content_as_str()
        self.lines = self.read_content_as_lines()
        self.atom_section = self.read_atoms_section_as_str()

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
            if "ATOM" or "HETATM" in line[0:6]:
                atoms_sect.append(line)
        return atoms_sect

    def read_atoms_section_as_str(self):
        return "\n".join(self.read_atoms_section())


def get_resname_from_line(line):
    return line[17:20]


def get_chain_from_line(line):
    return line[21:22]


def get_resnum_from_line(line):
    return line[22:26]


def set_resname_to_line(line, value):
    line[17:20] = value
    return line


def set_chain_to_line(line, value):
    line[21:22] = value
    return line


def set_resnum_to_line(line, value):
    line[22:26] = value
    return line




