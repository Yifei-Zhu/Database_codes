import h5py
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import defaultdict

from collections import Counter


def read_smiles_and_group_by_molecular_formula(hdf5_file):
    formula_to_smiles = defaultdict(list)

    with h5py.File(hdf5_file, 'r') as file:
        for dataset_name in file:
            dataset = file[dataset_name]
            smiles = dataset.attrs['SMILES']
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                formula = rdMolDescriptors.CalcMolFormula(mol)
                formula_to_smiles[formula].append(dataset_name)

    return formula_to_smiles

def group_by_molecular_formula(hdf5_file):
    formula_dict = defaultdict(list)
    with h5py.File(hdf5_file, 'r') as file:
        for dataset_name in file:
            dataset = file[dataset_name]
            labels = dataset['ground_state']['labels'][()][0].tolist()
            atomic_symbols = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F'}
            # 计算列表中每个元素的出现次数
            atom_counts = Counter(labels)
            # 生成分子式
            formula = ''
            for atomic_number in sorted(atom_counts.keys(), key=lambda x: (x != 6, x != 1, x)):
                symbol = atomic_symbols[atomic_number]
                count = atom_counts[atomic_number]
                formula += symbol + (str(count) if count > 1 else '')
            formula_dict[formula].append(dataset_name)
    return formula_dict

def test():
    hdf5_file = '/data/home/zhuyf/dataset_work/database/qm9_database/qm9_molecules.hdf5'
    formula_to_smiles = read_smiles_and_group_by_molecular_formula(hdf5_file)

    for formula, smiles_list in formula_to_smiles.items():
        print(f"{formula} : {smiles_list}")

if __name__ == '__main__':
    test()
