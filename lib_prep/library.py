import os
import glob
import subprocess
import multiprocessing as mp
from conf_variables import SCH_PATH

CMD = '{}/utilities/structconvert'.format(SCH_PATH) + " -isd {} -opdb {}"


class Library:
    def __init__(self, input_file):
        self.input_file = input_file
        self.content = self.read_content()
        self.filename, self.format = os.path.splitext(self.input_file)
        self.molecules = self.get_molecules()

    def read_content(self):
        with open(self.input_file) as infile:
            content = infile.read()
        return content

    def get_molecules(self):
        infile = self.content
        if self.format == ".sdf":
            list_of_mol = get_molecules_as_list(infile)
            molecules = {}
            for mol in list_of_mol:
                name = mol.split("\n")[0]
                counter = 1
                if name in molecules.keys():
                    print("Repeated names... Adding numerations")
                    while True:
                        new_name = name + "_{}".format(counter)
                        if new_name not in molecules.keys():
                            print("{} has been assigned to {}".format(name, new_name))
                            molecules[new_name] = mol
                            break
                        else:
                            counter += 1
                else:
                    molecules[name] = mol
            return molecules
        else:
            raise FileNotFoundError("Wrong format. Ensure that this an SDF file!")

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
        pool = mp.Pool(os.cpu_count())
        n_iter = len(all_sdfs)
        chunk_size = int((n_iter/200)/os.cpu_count())
        print("Transforming {} SDFs using {} cpus with a chunksize of {}...".format(n_iter, os.cpu_count(), chunk_size))
        pool.map(transform_sdf_to_pdb, all_sdfs, chunksize=chunk_size)


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


def transform_sdf_to_pdb(sdf):
    out_path = "/".join(sdf.split("/")[0:-2])+"/pdbs/"+os.path.splitext(os.path.basename(sdf))[0]+".pdb"
    command = CMD.format(sdf, out_path)
    subprocess.call(command.split())
