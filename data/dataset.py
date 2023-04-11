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
    with open(f'{zinc_code}.xyz', 'w') as xyz_file:
        xyz_file.write(f'{len(atoms)}\n')
        xyz_file.write(f'{zinc_code}\n')
        for atom in atoms:
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            xyz_file.write(f'{atom.GetSymbol()} {pos.x} {pos.y} {pos.z}\n')

def create_orca_input(zinc_code):
    with open(f'{zinc_code}.inp', 'w') as orca_input:
        orca_input.write(f'! RHF def2-SVP TightSCF PAL4\n')
        orca_input.write(f'*xyz 0 1\n')
        with open(f'{zinc_code}.xyz', 'r') as xyz_file:
            for line in xyz_file.readlines()[2:]:
                orca_input.write(line)
        orca_input.write(f'*\n')

def run_orca(zinc_code, orca_path='/Applications/orca/orca'):
    subprocess.run([orca_path, f'{zinc_code}.inp', f'-o', f'{zinc_code}.out'])

def extract_scf_time(zinc_code):
    with open(f'{zinc_code}.out', 'r') as orca_output:
        content = orca_output.read()
        scf_time_match = re.search(r'TOTAL SCF TIME\s+:\s+([\d.]+)s', content)
        if scf_time_match:
            scf_time = float(scf_time_match.group(1))
            return scf_time
        else:
            return None

def main():
    input_csv = 'substances.csv'
    output_csv = 'scf_times.csv'

    with open(input_csv, 'r') as csvfile:
        reader = pd.read_csv(csvfile)
        print(reader)
        scf_times = []

        for pos, mol in enumerate(reader["smiles"]):
            smiles = mol
            zinc_code = reader['zinc_id'][pos]
            print(f'Processing {zinc_code}...')
            try:
                smiles_to_xyz(smiles, zinc_code)
                create_orca_input(zinc_code)
                run_orca(zinc_code)

                scf_time = extract_scf_time(zinc_code)
                scf_times.append([zinc_code, scf_time])
            except:
                pass

    with open(output_csv, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ZINC_Code', 'SCF_Time'])
        writer.writerows(scf_times)

if __name__ == "__main__":
    main()
