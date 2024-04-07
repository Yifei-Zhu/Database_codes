import pandas as pd
from check1_OptedInchi import main_import_xyz as generateXYZ
from check1_OptedInchi import main_import_inchi as checkInchi
from check2_OptedRedundantInternalCoords import compare_ric_single as compareRicSingle
from check1_OptedInchi import xyz_to_inchi, xyz_to_smiles
from collections import defaultdict
import os
from check3_GraphBasedCheckDuplicate import main_import as checkDuplicateWithGraph

def _InchiBasedCheck(input_df, opted_index, new_xyz_path):
    different_file = f'{gau_path}/{label}_error_inchi.csv'
    opted_index['same_Inchi'], different_index = checkInchi(label, main_path, gau_path,new_xyz_path, input_df, different_file)
    return opted_index, different_index

def _RICBasedCheck(opted_index, index_list, new_xyz_path):
    df1 = pd.DataFrame(columns=['Index', 'Smiles_relaxed','InChI_relaxed'])
    nu_new=0
    for ind2 in index_list:
        gau_file = f'{gau_path}/{ind2}/{ind2}.log'
        s = compareRicSingle(gau_file)
        if s == True:
            opted_index['same_ric'].append(ind2)
        else:
            nu_new += 1
            xyz_file = f'{new_xyz_path}/{ind2}.xyz'
            inchi = xyz_to_inchi(xyz_file)
            smi = xyz_to_smiles(xyz_file)
            opted_index['new_strus'].append(ind2)
            df1 = df1.append({'Index': nu_new, 'Smiles_relaxedB': smi, 'InChI_relaxed': inchi}, ignore_index=True)

    base, extension = os.path.splitext(input_file)
    df1.to_csv(f'{base}_new{extension}')
    return opted_index

def getOptedGeoms(input_df):
    opted_index = defaultdict(list)
    new_xyz_path = f'{gau_path}/xyz_file_all'

    generateXYZ(gau_path, new_xyz_path, input_df)

    opted_index, different_index = _InchiBasedCheck(input_df, opted_index, new_xyz_path)
    # different_index = input_df

    opted_index = _RICBasedCheck(opted_index, different_index, new_xyz_path)
    return opted_index

def MAIN(input_df):
    opted_index = getOptedGeoms(input_df)
    combine_opted_index = opted_index['same_Inchi'] + opted_index['new_strus'] + opted_index['same_ric']
    combine_opted_index = checkDuplicateWithGraph(combine_opted_index, gau_path)
    with open(f'{main_path}/final_index.dat', 'w') as file:
        for item in combine_opted_index:
            file.write(f"{item}\n")


def getErrorList(gau_path, label):
    error_list_file = f'{gau_path}/error_list_{label}.dat'
    error_list = []
    with open(error_list_file, 'r') as f1:
        for line in f1:
            error_list.append(line.split(':'))
    return error_list


if __name__ == '__main__':
    global main_path, label, gau_path,input_file

    label = 'aromatic'
    # label = 'bond_but_nonaromatic'
    main_path = f'/data/home/zhuyf/dataset_work/database/qm9_database'
    gau_path = f'{main_path}/Gaussian'
    input_file = f'{main_path}/{label}.csv'
    error_list = getErrorList(gau_path, label)

    input_df = pd.read_csv(input_file)
    input_df = input_df[~input_df['Index'].isin(error_list)]['Index'][:100]

    MAIN(input_df)