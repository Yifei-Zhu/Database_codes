import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def remove_terminal_groups(mol):
    # 创建可编辑的分子
    mol_rw = Chem.RWMol(mol)
    # 找到端基原子
    terminal_atoms = [atom for atom in mol_rw.GetAtoms() if len(atom.GetNeighbors()) == 1]
    # 删除端基原子及其连接
    for atom in terminal_atoms:
        mol_rw.RemoveAtom(atom.GetIdx())
    # 将可编辑的分子转回普通分子
    mol_no_terminal = mol_rw.GetMol()

    # smiles_without_terminal = Chem.MolToSmiles(mol_no_terminal, isomericSmiles=True)
    return mol_no_terminal

def has_heteroatoms_without_terminal(mol):
    mol=remove_terminal_groups(mol)
    atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    return any(symbol not in ['C', 'H'] for symbol in atom_symbols)

def has_heteroatoms_in_rings(mol):
    ring_info = mol.GetRingInfo()
    for ring_atoms in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(atom).GetAtomicNum() not in [6, 1] for atom in ring_atoms):
            return True
    return False
        # return any(mol.GetAtomWithIdx(atom).GetAtomicNum() not in [6, 1] for atom in ring_atoms)

def is_fused_rings(mol):
    ring_info = mol.GetRingInfo()
    return len(ring_info.AtomRings()) >1 and any(len(set(ring1).intersection(set(ring2))) > 0 for ring1 in ring_info.AtomRings() for ring2 in ring_info.AtomRings())


def classify_molecule(mol):
    if mol:
        ring_info = mol.GetRingInfo()
        has_rings = ring_info.NumRings() > 0
        is_aromatic = Descriptors.NumAromaticRings(mol) > 0
        is_hetero_ring= has_heteroatoms_in_rings(mol)
        # is_hetero=any(symbol not in ['C', 'H'] for symbol in atom_symbols)
        is_hetero = has_heteroatoms_without_terminal(mol)
        is_fused = is_fused_rings(mol)

        if has_rings:
            if is_aromatic and is_hetero_ring:
                return 1
            elif is_aromatic and not is_hetero_ring:
                return 2
            elif is_fused and is_hetero_ring:
                return 3
            elif is_fused and not is_hetero_ring:
                return 4
            elif is_hetero_ring:
                return 5
            else:
                return 6
        else:
            if is_hetero:
                return 8
            else:
                return 7
    else:
        return 9

def classify_main(df, col_label1, col_label2):
    compounds_type = {
        1: 'heteroaromatics',
        2: 'aromatics',
        3: 'fused heterocycles',
        4: 'fused carbocycles',
        5: 'heterocycles',
        6: 'carbocycles',
        7: 'carboacyclic',
        8: 'heteroacyclic',
        9: 'error'
    }

    def classify_row(row):
        mol = None
        if row[col_label1] != '1':
            mol = Chem.MolFromSmiles(row[col_label1])

        if mol is None:
            mol = Chem.MolFromInchi(row[col_label2])

        if mol is not None:
            return compounds_type[classify_molecule(mol)]
        else:
            return None

    df['CompoundType'] = df.apply(classify_row, axis=1)
    return df

def test():
    compounds_type = {1:'heteroaromatics',2:'aromatics',3:'fused heterocycles',4:'fused carbocycles',5:'heterocycles',6:'carbocycles',7:'carboacyclic',8:'heteroacyclic'}
    smiles_list = ['c1ccncc1', 'c1ccccc1', 'C1CCOC1', 'C1CCC2CCCCC2C1', 'C1CCOC1', 'C1CCC(C1)O', 'CCOC(=O)C']
    print(compounds_type[classify_molecule('C(CO)O')])

if __name__=='__main__':
    test()

    # code_path = '/data/home/zhuyf/dataset_work/database/analysis'
    # main_path = '/data/home/zhuyf/dataset_work/database/qm9_database'
    # col_label = 'Smiles_GDB'

    # smi_file = f'{main_path}/double_and_triple.csv'
    # smi_file = f'{main_path}/a.csv'
    # df = pd.read_csv(smi_file, index_col=0)

    # classify_main(df)

