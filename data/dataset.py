import csv
import os
import subprocess
import re
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

def smiles_to_xyz(smiles, zinc_code):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    atoms = mol.GetAtoms()
    with open(f'./{zinc_code}/{zinc_code}.xyz', 'w') as xyz_file:
        xyz_file.write(f'{len(atoms)}\n')
        xyz_file.write(f'{zinc_code}\n')
        for atom in atoms:
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            xyz_file.write(f'{atom.GetSymbol()} {pos.x} {pos.y} {pos.z}\n')

def create_orca_input(zinc_code):
    with open(f'./{zinc_code}/{zinc_code}.inp', 'w') as orca_input:
        orca_input.write(f'!RHF def2-SVP TightSCF PAL16 \n')
        orca_input.write(f'*xyz 0 1\n')
        with open(f'./{zinc_code}/{zinc_code}.xyz', 'r') as xyz_file:
            for line in xyz_file.readlines()[2:]:
                orca_input.write(line)
        orca_input.write(f'*\n')

def run_orca(zinc_code, orca_path='/opt/orca-5.0.3/orca'):
    command = f"{orca_path} ./{zinc_code}/{zinc_code}.inp | tee ./{zinc_code}/{zinc_code}.out"
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True, text=True)

def read_output_file(filename):
    with open(filename, 'r') as file:
        content = file.read()
    return content

def extract_total_time(zinc_code):
    # Define a regular expression pattern to search for "Total time"
    content = read_output_file(f"./{zinc_code}/{zinc_code}.out")

    pattern = r"Total time\s*\.+\s*([\d.]+) sec"

    # Search for the pattern in the content
    match = re.search(pattern, content)

    if match:
        total_time = float(match.group(1))
        return total_time
    else:
        return None

def main():
    input_csv = 'substances.csv'
    output_csv = open('scf_times.csv', "a")
    output_csv.write("ZINC_ID,SMILES,SCF\n")
    with open(input_csv, 'r') as csvfile:
        reader = pd.read_csv(csvfile)
        print(reader)
        scf_times = []

        for pos, mol in enumerate(reader["smiles"]):
            smiles = mol
            zinc_code = reader['zinc_id'][pos]
            print(f'Processing {zinc_code}...')

            os.makedirs(f"./{zinc_code}", exist_ok=True)
            try:
                smiles_to_xyz(smiles, zinc_code)
                create_orca_input(zinc_code)
                run_orca(zinc_code)
                print("DONE!")

                scf_time = extract_total_time(zinc_code)
                print(scf_time)
                scf_times.append([zinc_code, scf_time])

                output_csv.write(f"{zinc_code},{smiles},{scf_time}\n")
                output_csv.flush()
            except:
                print("FAILED!")
                pass

    input_csv.close()
    output_csv.close()
if __name__ == "__main__":
    main()
