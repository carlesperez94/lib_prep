import os
import glob
import subprocess
import multiprocessing as mp
from rdkit import Chem
from lib_prep.conf_variables import SCH_PATH

CMD = '{}/utilities/structconvert'.format(SCH_PATH) + " -isd {} -opdb {}"


class LibrarySDF:
    def __init__(self, input_file):
        self.input_file = input_file
        self.content = self.read_content()
        self.filename, self.format = os.path.splitext(self.input_file)
        self.molecules = {}
        self.outpath = None

    def read_content(self):
        with open(self.input_file) as infile:
            content = infile.read()
        return content

    def get_molecules(self):
        infile = self.content
        if self.format == ".sdf":
            list_of_mol = get_molecules_as_list(infile)
            for mol in list_of_mol:
                name = mol.split("\n")[0]
                counter = 1
                if name in self.molecules.keys():
                    print("Repeated names... Adding numerations")
                    while True:
                        new_name = name + "_{}".format(counter)
                        if new_name not in self.molecules.keys():
                            print("{} has been assigned to {}".format(name, new_name))
                            self.molecules[new_name] = mol
                            break
                        else:
                            counter += 1
                else:
                    self.molecules[name] = mol
            return self.molecules
        else:
            raise FileNotFoundError("Wrong format. Ensure that this an SDF file!")

    def get_molecules_as_rdkit(self, remove_h=False):
        molecules = Chem.SDMolSupplier(self.input_file, removeHs=remove_h)
        for mol in molecules:
            yield mol

    def export_sdf_to_sdfs(self, out_path):
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        sdfs_out_path = os.path.join(out_path, "sdfs/")
        if not os.path.exists(sdfs_out_path):
            os.mkdir(os.path.join(sdfs_out_path))
        for name, molecule in self.molecules.items():
            with open(os.path.join(out_path, sdfs_out_path+name+".sdf"), "w") as out_file:
                out_file.write(molecule)

    def export_sdf_to_pdbs(self, out_path):
        self.export_sdf_to_sdfs(out_path)
        all_sdfs = glob.glob("{}/sdfs/*.sdf".format(out_path))
        pdb_path = out_path+"/pdbs"
        os.mkdir(pdb_path)
        configure_multiprocessing(all_sdfs, transform_sdf_to_pdb)

    def set_outpath(self, out_path):
        self.outpath = out_path

    def get_sdf_filenames(self):
        molecules_names = []
        molecules_split = self.content.split("$$$$\n")
        name_1 = molecules_split[0].split("\n")[0]
        molecules_names.append(name_1)
        for molecule in molecules_split[1:-1]:
            molecule_lines = molecule.split("\n")
            molname = molecule_lines[0]
            molecules_names.append(molname)
        return molecules_names


    def export_sdf_to_pdb_rdkit(self, out_path, molname_property="Molecule Name", get_sdf_name=True):
        self.set_outpath(out_path)
        if not os.path.exists(self.outpath):
            os.mkdir(self.outpath)
        pdb_path = out_path + "/pdbs"
        if not os.path.exists(pdb_path):
            os.mkdir(pdb_path)
        molecule_names = {}
        molecule_counter = 1
        if get_sdf_name:
            mol_names = self.get_sdf_filenames()
        for n, mol in enumerate(self.get_molecules_as_rdkit()):
            if not get_sdf_name:
                try:
                    mol_name = mol.GetPropsAsDict()[molname_property]
                except KeyError:
                    mol_name = "MOL{:04d}".format(molecule_counter)
                    molecule_counter += 1
            else:
                mol_name = mol_names[n]
            add_to_dictionary_and_count(dictionary=molecule_names, entrance=mol_name)
            new_name = "{}_{}.pdb".format(mol_name, molecule_names[mol_name])
            new_filename = os.path.join(pdb_path, new_name)
            writer = Chem.PDBWriter(new_filename)
            writer.write(mol)
            print("{} Succesfully transformed!".format(new_filename))


def find_molecules_limit_sdf(file_content):
    lines = file_content.split("\n")
    molecules_limit = []
    for n in range(1000000000):
        try:
            if "$$$$" in lines[n]:
                molecules_limit.append(n+1)
        except IndexError:
            print("Analysis finished")
            return molecules_limit


def get_molecules_as_list(file_content):
    limits = find_molecules_limit_sdf(file_content)
    lines = file_content.split("\n")
    molecules = []
    for n in range(0, len(limits)):
        if n == 0:
            pass
        else:
            mol = lines[limits[n-1]:limits[n]]
            molecules.append("\n".join(mol))
    return molecules


def add_to_dictionary_and_count(dictionary, entrance):
    if entrance in dictionary:
        dictionary[entrance] = dictionary[entrance] + 1
    else:
        dictionary[entrance] = 0


def transform_sdf_to_pdb(sdf):
    out_path = "/".join(sdf.split("/")[0:-2])+"/pdbs/"+os.path.splitext(os.path.basename(sdf))[0]+".pdb"
    command = CMD.format(sdf, out_path)
    subprocess.call(command.split())


def configure_multiprocessing(input_list_task, function_to_apply):
    pool = mp.Pool(os.cpu_count())
    n_iter = len(input_list_task)
    chunk_size = int((n_iter / 200) / os.cpu_count())
    print("{} tasks using {} cpus with a chunksize of {}...".format(n_iter, os.cpu_count(), chunk_size))
    pool.map(function_to_apply, input_list_task, chunksize=chunk_size)
