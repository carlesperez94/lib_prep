import lib_prep.convert_sdfs2pdb as s2p
import lib_prep.pdb_modifier as pdbm
import lib_prep.FragmentTools.prepare_to_frag as p2f
import glob
import os
import shutil

RESULT_TEST = ["results_testing_convert/pdbs/F6513-2497_0.pdb", "results_testing_convert/pdbs/F6513-2497_1.pdb",
               "results_testing_convert/pdbs/F6513-2497_2.pdb",  "results_testing_convert/pdbs/F6540-1156_0.pdb",
               "results_testing_convert/pdbs/F6540-1156_1.pdb"]


class Convert2SDFTest:

    CHAIN = "Z"
    RESNAME = "TST"
    RESNUM = 999

    def test_files_creation(self):
        s2p.main(sdf_input="libraries/testing_set.sdf", output_folder="results_testing_convert",
                 chain=self.CHAIN, resname=self.RESNAME, resnum=self.RESNUM)
        result_files = glob.glob("results_testing_convert/pdbs/*")
        assert sorted(result_files) == RESULT_TEST

    def test_file_content(self):
        result_files = glob.glob("results_testing_convert/pdbs/*")
        chains = []
        resnames = []
        resnums = []
        for pdb in result_files:
            pdb_class = pdbm.PDB(pdb)
            for line in pdb_class.read_atoms_section():
                chains.append(pdbm.get_chain_from_line(line))
                resnames.append(pdbm.get_resname_from_line(line).strip())
                resnums.append(pdbm.get_resnum_from_line(line).strip())
        assert list(set(chains)) == [self.CHAIN] and list(set(resnames)) == [self.RESNAME] \
               and list(set(resnums)) == [str(self.RESNUM)]


class Prepare2FragTest:
    SERIE_FILE = "serie_testing.conf"
    HV_ATOM = "C1"

    def test_file_creation(self):
        p2f.main(pdb_complex="None", heavy_atom_pdb_complex=self.HV_ATOM, lib_path="results_testing_convert/pdbs",
                 mode="first-occurrence", out_file=self.SERIE_FILE)
        assert os.path.exists(self.SERIE_FILE)

    def test_file_content(self):
        with open(self.SERIE_FILE) as ser_file:
            content_as_lines = ser_file.readlines()
            assert len(content_as_lines) == 5
            for line in content_as_lines:
                pdb_path = line.split()[0]
                pdb_relative_path = pdb_path.split("/")[-3:]
                pdb_relative_path = "/".join(pdb_relative_path)
                assert pdb_relative_path in RESULT_TEST


def testing_cleaning():
    if os.path.exists("results_testing_convert"):
        shutil.rmtree("results_testing_convert")
    if os.path.exists("serie_testing.conf"):
        os.remove("serie_testing.conf")