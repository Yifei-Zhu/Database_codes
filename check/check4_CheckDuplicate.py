import os
import pandas as pd
import numpy as np
import h5py
import json

from sub_check41_Classify import group_by_molecular_formula as classify_mols
from sub_check42_GraphCheck import main_import as graph_check
from sub_check42_Smiles_and_Inchi import check_by_smiles, check_by_inchi,get_duplicate_file_inchi,get_duplicate_file_smiles

def main_test():
    global main_path, label, gau_path,input_file

    label = 'double_and_triple'
    main_path = f'/data/home/zhuyf/dataset_work/database/qm9_database'
    gau_path = f'{main_path}/Gaussian'
    hdf5_file = '/data/home/zhuyf/dataset_work/database/qm9_database/qm9_molecules.hdf5'
    formula_dict = classify_mols(hdf5_file)
    duplicate = {}
    for formula, ind_list in formula_dict.items():
        duplicate[formula] = graph_check(ind_list, gau_path)
    print(duplicate)

def remove_duplicate_csv(final_csv, duplicate_list):
    df = pd.read_csv(final_csv, low_memory=False)
    df_filtered = df[~df['Index'].isin(duplicate_list)]
    df_filtered = df_filtered.loc[:, ~df_filtered.columns.str.contains('Unnamed')]
    df_filtered = df_filtered.reset_index()
    #df_filtered.index+=1
    df_filtered.to_csv(f"final_all_remove_duplicate_rev.csv")


def remove_duplicate_hdf5(dup_list, final_hdf5_file):
    with h5py.File(final_hdf5_file, 'a') as file:
        for dataset_name in dup_list:
            if dataset_name in file:
                del file[dataset_name]


def remove_duplicate(duplicate_list, final_csv, final_hdf5_file):
    # remove_duplicate_hdf5(duplicate_list, final_hdf5_file)
    remove_duplicate_csv(final_csv, duplicate_list)


def get_formula_dict(check_path, final_hdf5_file):
    formula_dict = classify_mols(final_hdf5_file)
    with open(f'{check_path}/formula_dict.json', 'w') as f:
        json.dump(formula_dict, f)


def check_duplicate_graph(check_path, datasets_info):
    with open(f'{check_path}/formula_dict.json', 'r') as f:
        formula_dict = json.load(f)
    duplicate_list = []
    for formula, ind_list in formula_dict.items():
        dup_list = graph_check(ind_list, check_path, datasets_info)
        if dup_list:
            duplicate_list.extend(dup_list)
        print(f'{formula}:{len(dup_list)}')
    # np.save(f'{check_path}/duplicate_all.npy', duplicate_list)
    return duplicate_list

def main():
    check_path = '/data/home/zhuyf/dataset_work/database/checkOptedGeoms'
    final_csv = f"{check_path}/final_all.csv"
    final_hdf5_file = f"{check_path}/final_all_with_duplicates.hdf5"
    datasets = np.load(f'{check_path}/datasets_dict.npy',allow_pickle=True).item()

    choose = str(input('1 - get_formula_dict \n2 - get_duplicate_list \n3 - remove_duplicate\n')).strip()

    if choose == '1' :
        get_formula_dict(check_path, final_hdf5_file)
    elif choose == '2':
        choose_i = int(input('1 - Graph-based \n2 - Smiles-based \n3 - Inchi-based\n'))
        if choose_i == 1:
            duplicate_list = check_duplicate_graph(check_path, datasets)
            json_file=f'{check_path}/duplicate_list_graph.json'
        elif choose_i ==2:
            duplicate_list = check_by_smiles(final_hdf5_file)
            json_file = f'{check_path}/duplicate_list_smiles_rev.json'
        elif choose_i == 3:
            duplicate_list = check_by_inchi(final_hdf5_file)
            json_file = f'{check_path}/duplicate_list_inchi_rev.json'

        with open(json_file, 'w') as jsonfile:
            json.dump(duplicate_list, jsonfile)

    elif choose == '3':
        os.chdir(check_path)
        json_file = f'{check_path}/duplicate_list_smiles_rev.json'
        with open(json_file, 'r') as jsonfile:
            tmp1 = json.load(jsonfile)
        json_file = f'{check_path}/duplicate_list_inchi_rev.json'
        with open(json_file, 'r') as jsonfile:
            tmp2 = json.load(jsonfile)
        dup_list = set(tmp1+tmp2)

        # json_file = f'{check_path}/duplicate_list_SMILESandINCHI.json'
        # with open(json_file, 'r') as jsonfile:
        #     dup_list = json.load(jsonfile)

        # a=set(dup_list2)^set(dup_list1)
        # print(len(a))
        remove_duplicate(dup_list, final_csv, final_hdf5_file)


if __name__ == '__main__':
    main()
