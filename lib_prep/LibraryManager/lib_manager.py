import glob
from lib_prep import pdb_modifier


class LibraryChecker:
    def __init__(self, path_to_library):
        self.path_to_library = path_to_library
        self.format = self.check_library_format()

    def get_files(self):
        library_files = glob.glob("{}/*".format(self.path_to_library))
        return library_files

    def check_library_format(self):
        pdb = True
        mae = True
        sdf = True
        for file in self.get_files():
            if not file.endswith(".pdb"):
                pdb = False
            if not file.endswith(".mae"):
                mae = False
            if not file.endswith(".sdf"):
                sdf = False

        if pdb == True and mae == False and sdf == False:
            return ".pdb"

        elif mae == True and pdb == False and sdf == False:
            return ".mae"

        elif sdf == True and pdb == False and mae is False:
            return ".sdf"
        else:
            print("PDB:{}, MAE:{}, SDF:{}".format(pdb, mae, sdf))
            raise IOError("Formats are mixed in the input folder. Please, ensure that all files have the same format.")

    def prepare_library(self, chain, resname, resnum):
        if self.format == ".pdb":
            for pdb in self.get_files():
                pdb_object = pdb_modifier.PDB(in_pdb=pdb, chain=chain, resname=resname, resnum=resnum)
                pdb_object.modify_pdb()
                pdb_object.overwrite_pdb()