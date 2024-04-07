from openbabel import openbabel as ob
from sub_code_get_GS_opt_props import main_xyz
import os
from rdkit import Chem
import pandas as pd
import csv
from openbabel import pybel


def generate_xyz_from_log(input_df):
    for index_i in input_df:
        print(index_i)
        file_opt = f'{gau_path}/{index_i}/{index_i}.log'
        xyz_file = f'{new_xyz_path}/{index_i}.xyz'
        main_xyz(file_opt, xyz_file, 'opt+freq')


def get_inchi_after_opt(input_df,new_xyz_path,new_inchi_file):
    df2 = pd.DataFrame(columns=['Index', 'InChI'])
    for index_i in input_df:
        print(index_i)
        xyz_file = f'{new_xyz_path}/{index_i}.xyz'
        inchi = xyz_to_inchi(xyz_file)
        # xyz_to_mol(xyz_file, mol_file)
        # inchi = get_new_inchi(mol_file)

        df2 = df2.append({"Index": index_i, "InChI": inchi}, ignore_index=True)
    df2.index+=1
    df2.to_csv(new_inchi_file)

def get_inchi_before_opt(input_df,old_xyz_path,old_inchi_file):
    df2 = pd.DataFrame(columns=['Index', 'InChI'])
    for index_i in input_df:
        print(index_i)
        xyz_file = f'{old_xyz_path}/{index_i}.xyz'
        inchi = xyz_to_inchi(xyz_file)
        # xyz_to_mol(xyz_file, mol_file)
        # inchi = get_new_inchi(mol_file)

        df2 = df2.append({"Index": index_i, "InChI": inchi}, ignore_index=True)
    df2.index+=1
    df2.to_csv(old_inchi_file)

# def xyz_to_inchi(input_xyz_file):
#     obConversion = ob.OBConversion()
#     obMolecule = ob.OBMol()
#     obConversion.SetInAndOutFormats("xyz", "inchi")
#     obConversion.ReadFile(obMolecule, input_xyz_file)
#     inchi = obConversion.WriteString(obMolecule)
#     return inchi

def xyz_to_inchi(input_xyz_file):
    # 使用pybel读取XYZ文件
    molecule = pybel.readfile("xyz", input_xyz_file).__next__()
    # 将分子转换为InChI字符串
    inchi = molecule.write("inchi")
    return inchi

def xyz_to_smiles(input_xyz_file):
    molecule = pybel.readfile("xyz", input_xyz_file).__next__()
    smi = molecule.write(format='smi').strip()
    return smi

def compare_inchi(old_inchi_file,new_inchi_file, different_file):

    inchi_old = pd.read_csv(old_inchi_file)
    inchi_new = pd.read_csv(new_inchi_file)

    merged_df = pd.merge(inchi_old, inchi_new, on='Index', suffixes=('_old', '_new'), how='inner')

    different_index = merged_df[merged_df['InChI_old'] != merged_df['InChI_new']]['Index']
    same_index = merged_df[merged_df['InChI_old'] == merged_df['InChI_new']]['Index']

    different_index.to_csv(different_file)

    return same_index.tolist(), different_index.tolist()

def main_import_xyz(gau_path, new_xyz_path, input_df):
    os.makedirs(new_xyz_path, exist_ok=True)
    generate_xyz_from_log(gau_path, new_xyz_path, input_df)

def main_import_inchi(label, main_path, gau_path,new_xyz_path, input_df, different_file):
    new_inchi_file = f'{gau_path}/{label}_new_inchi.csv'
    old_inchi_file = f'{main_path}/old_inchi.csv'
    get_inchi_before_opt(input_df,old_xyz_path,old_inchi_file)
    get_inchi_after_opt(input_df,new_xyz_path,new_inchi_file)
    same_index, different_index = compare_inchi(old_inchi_file,new_inchi_file,different_file)
    return same_index,different_index


if __name__ == '__main__':
    global main_path, label, gau_path

    #label = 'double_and_triple'
    #main_path = f'/data/home/zhuyf/dataset_work/database/qm9_database'
    #gau_path = f'{main_path}/Gaussian'
    #input_file = f'{main_path}/{label}.csv'
    #error_list_file = f'{gau_path}/error_list_{label}.csv'
    #error_list = pd.read_csv(error_list_file)['Index']
    #input_df = pd.read_csv(input_file)
    #input_df = input_df[~input_df['Index'].isin(error_list)]['Index']

    label=''
    #label='excited'
    main_path = '/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation'
    prefix = 'gdb11_'
    gau_path = f'{main_path}/calculation_{label}/{prefix}Gaussian_{label}'
    input_file = f'{main_path}/new_mbkmeans_final_sample_1w.csv'
    error_list_file = f'{gau_path}/error_list_{label}.csv'
    error_list = pd.read_csv(error_list_file)['Index']
    input_df = pd.read_csv(input_file)
    input_df = input_df[~input_df['total_index_from_random_sample'].isin(error_list)]['total_index_from_random_sample']
    old_xyz_path = '/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation/calculation_/gdb11_initial_xyz_'


    choose = int(input('0 - generate xyz_file from log file\n1 - InChI-based check\n'))
    global new_xyz_path
    new_xyz_path = f'{gau_path}/xyz_file_ala'

    if choose == 0:
        os.makedirs(new_xyz_path, exist_ok=True)
        generate_xyz_from_log(input_df)

    elif choose == 1:
        new_inchi_file = f'{gau_path}/{label}_new_inchi.csv'
        old_inchi_file = f'{main_path}/old_inchi.csv'
        different_file = f'{gau_path}/{label}_error_inchi.csv'

        get_inchi_before_opt(input_df,old_xyz_path,old_inchi_file)
        get_inchi_after_opt(input_df,new_xyz_path,new_inchi_file)
        compare_inchi(old_inchi_file,new_inchi_file,different_file)
