import os
import re
import sys
import numpy as np
import pandas as pd
import json
import pickle
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from collections import defaultdict
from openbabel import pybel
import periodictable


def xyz_to_smiles(input_xyz_file):
    molecule = pybel.readfile("xyz", input_xyz_file).__next__()
    smiles = molecule.write("smiles", opt={"h": None})
    return smiles.split()[0]

def xyz_to_inchi(input_xyz_file):
    molecule = pybel.readfile("xyz", input_xyz_file).__next__()
    inchi = molecule.write("inchi")
    return inchi
'''
def generate_inchi_from_xyz(xyz_coordinates, atomic_symbols):
    # 创建 RDKit 分子对象

    mol = Chem.MolFromSmiles('')  # 创建一个空分子对象
    mol = Chem.RWMol(mol)

    # 添加原子和坐标信息到分子对象
    for symbol in atomic_symbols:
        atom = Chem.Atom(symbol)
        mol.AddAtom(atom)

    conf = Chem.Conformer(len(xyz_coordinates))
    for i, coord in enumerate(xyz_coordinates):
        x, y, z = map(float, coord)
        conf.SetAtomPosition(i, (x, y, z))
    conf = Chem.Conformer(len(atomic_symbols))
    mol.AddConformer(conf)

    for i in range(len(atomic_symbols)):
        for j in range(i+1, len(atomic_symbols)):
            bond_length = mol.GetConformer().GetAtomPosition(i).Distance(mol.GetConformer().GetAtomPosition(j))
            mol.AddBond(i, j, Chem.BondType.SINGLE)
            mol.GetBondBetweenAtoms(i, j).SetProp('Length', str(bond_length))

    smi = Chem.MolToSmiles(mol)
    inchi = Chem.MolToInchi(mol)

    return smi, inchi
'''
# def remove_isolated_hydrogens(molecule):
#     # 查找孤立的氢原子
#     isolated_hydrogens = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 1 and atom.GetDegree() == 0]
#     mol = Chem.RWMol(molecule)
#     # 逆序删除孤立氢原子
#     for atom_idx in reversed(isolated_hydrogens):
#         mol.RemoveAtom(atom_idx)
#     return mol


# def xyz_to_smi_and_inchi(xyz_file):
#     atomic_labels, coords = readXYZ_single(xyz_file)
#     atomic_numbers = list(map(symbol_to_atomic_number_rdkit, atomic_labels))
#     inchi = generate_smi_and_inchi_from_xyz(coords, atomic_numbers)
#     return inchi

def save_smiles_inchi_to_file(smiles, inchi, output_file):
    with open(output_file, 'w') as file:
        file.write(f'SMILES:\n {smiles}\n')
        file.write(f'InChI:\n {inchi}\n')

def XyzToStrs_single(ind):
    print(ind)
    path_i = f'{gau_path}/{ind}'
    xyz_file = f'{path_i}/{ind}.xyz'
    SmiAndInchi = f'{path_i}/SmiAndInchi.dat'
    elements, coords = readXYZ_single(xyz_file)
    smi = xyz_to_smiles(xyz_file)
    inchi = xyz_to_inchi(xyz_file)
    # smi, inchi = generate_inchi_from_xyz(coords,elements)
    save_smiles_inchi_to_file(smi, inchi, SmiAndInchi)


def mainXyzToStrs(source_csv):
    df_ind = pd.read_csv(source_csv)['Index']
    for ind in df_ind:
        if ind != 129152:
            XyzToStrs_single(ind)


def save_bin(path_i, bin_data):
    file_name = f'{path_i}/ecfp_data.bin'
    with open(file_name, 'wb') as file:
        pickle.dump(bin_data, file)

def saveNpy(path_i, data):
    file_name = f'{path_i}/ecfp_data.npy'
    np.save(file_name, data)
    
def cal_ECFP(mol):
    ecfp4 = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    ecfp4_string = ecfp4.ToBitString()
    return ecfp4_string

def generateECFP_single(path_i, smi,inchi):
    try:
        mol = Chem.MolFromSmiles(smi)
        ecfp = cal_ECFP(mol)
    except:
        mol = Chem.MolFromInchi(inchi)
        ecfp = cal_ECFP(mol)
    saveNpy(path_i, ecfp)

def main_generaterECFP(ind):
    path_i = f'{gau_path}/{ind}'
    smi_file = f'{path_i}/SmiAndInchi.dat'
    with open(smi_file, 'r') as f1:
        lines=f1.readlines()
        smi = lines[1].split()[0].strip()
        inchi = lines[4].strip()
    generateECFP_single(path_i, smi, inchi)

def generateCM_single(atomic_labels, coords):
    atomic_numbers = list(map(symbol_to_atomic_number_rdkit, atomic_labels))
    num_atoms=len(atomic_labels)
    coulomb_matrix = np.zeros((num_atoms, num_atoms))
    for i in range(num_atoms):
        for j in range(num_atoms):
            if i == j:
                coulomb_matrix[i, j] = 0.5 * atomic_numbers[i]**2.4
            else:
                Zi, Zj = atomic_numbers[i], atomic_numbers[j]
                rij = np.linalg.norm(coords[i] - coords[j])
                coulomb_matrix[i, j] = Zi * Zj / rij

    return coulomb_matrix

def symbol_to_atomic_number_rdkit(element_symbol):
    element = periodictable.elements.symbol(element_symbol)
    atomic_number = element.number
    return atomic_number

def readXYZ_single(xyz_file):
    elements=[]
    with open(xyz_file, 'r') as f1:
        num_atoms = int(f1.readline().strip())
        f1.readline()
        coords=np.zeros((num_atoms, 3))
        for i in range(num_atoms):
            line = f1.readline().split()
            element = line[0]
            x, y, z = map(float, line[1:4])
            elements.append(element)
            coords[i] = [x, y, z]
    return elements, coords

def main_generateCM(ind):
    path_i = f'{gau_path}/{ind}'
    xyz_file = f'{path_i}/{ind}.xyz'
    atomic_labels, coords = readXYZ_single(xyz_file)
    coulomb_matrix=generateCM_single(atomic_labels, coords)
    np.savetxt(f'{path_i}/cm.dat',coulomb_matrix)


def main_generaterDescriptors(source_csv):
    df_ind = pd.read_csv(source_csv)['Index']
    # df_fps = pd.DataFrame(columns=['Index', 'Smiles_GDB','MACCS', 'ECFP'])
    for ind in df_ind:
        if ind != 129152:
            # XyzToStrs_single(ind)
            main_generaterECFP(ind)
            # main_generateCM(ind)

def main():
    global code_path, main_path, gau_path

    code_path = '/data/home/zhuyf/dataset_work/database/analysis'
    main_path = '/data/home/zhuyf/dataset_work/database/qm9_database'
    gau_path = f'{main_path}/Gaussian'

    smi_file = f'{main_path}/double_and_triple.csv'
    # smi_file = f'{main_path}/a.csv'

    # mainXyzToStrs(smi_file)
    main_generaterDescriptors(smi_file)

def main_test_smi_inchi():
    global code_path, main_path, gau_path

    code_path = '/data/home/zhuyf/dataset_work/database/analysis'
    main_path = '/data/home/zhuyf/dataset_work/database/qm9_database'
    gau_path='/data/home/zhuyf/dataset_work/database/qm9_database/database/xyz_file_all'
    smi_file = f'{main_path}/double_and_triple.csv'
    # smi_file = f'{main_path}/a.csv'
    df_ind = pd.read_csv(smi_file)
    for ind in df_ind['Index']:
        path_i = f'{gau_path}'
        xyz_file = f'{path_i}/{ind}.xyz'
        try:
            smi = xyz_to_smiles(xyz_file)
            inchi = xyz_to_inchi(xyz_file)
        except:
            df_ind.loc[df_ind['Index'] == ind, 'Smiles_relaxed'] = smi.split()[0]
            df_ind.loc[df_ind['Index'] == ind, 'InChI_relaxed'] = inchi.replace("\n", "")
    df_ind.to_csv(f'{code_path}/test.csv')

if __name__ == '__main__':
    main()
    # main_test_smi_inchi()
    print('Done!\n')
