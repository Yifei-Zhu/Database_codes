import h5py
import pandas as pd
import numpy as np
import tempfile
import subprocess as sub
import os

from openbabel import pybel
from rdkit import Chem

from sub_code_get_GS_opt_props import main_xyz as get_GS_opt_xyz

def generate_final_csv_each(target):
    target_name = target[0]
    target_dict = target[1]
    df1 = pd.read_csv(f"{target_dict['main_path']}/{target_dict['csv_file']}")
    df2 = pd.read_csv(f"{target_dict['gs_dir']}/{target_dict['gs_error']}")
    df3 = pd.read_csv(f"{target_dict['es_dir']}/{target_dict['es_error']}")

    primary_index='total_index_from_random_sample'
    secondary_index='Index'
    index_to_use = primary_index if primary_index in df1.columns else secondary_index
    filtered_df = filter_dataframe(df1, df2, df3, index_to_use)

    prefix_df = add_prefix(target_dict, filtered_df, index_to_use)
    if "Unnamed: 0" in prefix_df.columns:
        prefix_df.drop("Unnamed: 0", axis=1, inplace=True)
    prefix_df.reset_index(drop=True, inplace=True)
    prefix_df.index += 1
    prefix_df.to_csv(f"{target_dict['main_path']}/final_{target_name}.csv", index=False)


def filter_dataframe(df1, df2, df3, index):
    filtered_df = df1[~df1[index].isin(df2['Index'])]
    filtered_df = filtered_df[~filtered_df[index].isin(df3['Index'])]
    return filtered_df

def add_prefix(target_dict, df, index):
    df['Index'] = str(target_dict['prefix']) + df[index].astype(str)
    return df


def merge_final_csv(datasets, final_csv):
    columns = ['Index', 'Smiles', 'InChI', 'HeavyAtomCount', 'RingNumber', 'CompoundType']
    df_combined = pd.DataFrame(columns=columns)
    for target_name, target_dict in datasets.items():
        df = pd.read_csv(f"{target_dict['main_path']}/final_{target_name}.csv")
        df_filtered = df.loc[:, df.columns.intersection(columns)]  # 筛选出需要的列
        df_combined = pd.concat([df_combined, df_filtered], ignore_index=True)

    # 确保最终的DataFrame列的顺序正确
    df_combined = df_combined[columns]
    df_combined.to_csv(final_csv)
    return df_combined

def replace_smiles_inchi_single(prefix, hdf5_file, df):
    # Initialize empty lists to store updates
    updates = {
        'Index': [],
        'Smiles_pybel': [],
        'InchI_pybel': [],
        'Smiles_rdkit': [],
        'InchI_rdkit': []
    }

    with h5py.File(hdf5_file, 'r') as file:
        for dataset_name in file:
            dataset = file[dataset_name]
            prefixed_name = prefix + dataset_name
            if dataset.attrs['SMILES_PYBEL'] == '1' and dataset.attrs['INCHI_PYBEL'] == '1' and dataset.attrs['SMILES_RDKIT'] == '1' and dataset.attrs['INCHI_RDKIT'] == '1':
                continue
            # Append data to lists
            updates['Index'].append(prefixed_name)
            updates['Smiles_pybel'].append(dataset.attrs['SMILES_PYBEL'])
            updates['InchI_pybel'].append(dataset.attrs['INCHI_PYBEL'])
            updates['Smiles_rdkit'].append(dataset.attrs['SMILES_RDKIT'])
            updates['InchI_rdkit'].append(dataset.attrs['INCHI_RDKIT'])

    # Create a temporary DataFrame from the updates
    temp_df = pd.DataFrame(updates).set_index('Index')

    # Update the original DataFrame with the temporary one
    for column in temp_df.columns:
        df.update(temp_df[column], overwrite=True)

    return df

def replace_smiles_inchi_total(hdf5_file, df):
    with h5py.File(hdf5_file, 'r') as file:
        for dataset_name in file:
            print(dataset_name)
            dataset = file[dataset_name]
            smi1 = dataset.attrs['SMILES_PYBEL']
            inchi1 = dataset.attrs['INCHI_PYBEL']
            smi2 = dataset.attrs['SMILES_RDKIT']
            inchi2 = dataset.attrs['INCHI_RDKIT']
            df.loc[df['Index'] == dataset_name, 'Smiles_pybel'] = smi1
            df.loc[df['Index'] == dataset_name, 'InchI_pybel'] = inchi1
            df.loc[df['Index'] == dataset_name, 'Smiles_rdkit'] = smi2
            df.loc[df['Index'] == dataset_name, 'InchI_rdkit'] = inchi2
    return df


def replace_smiles_inchi(datasets, final_csv, final_hdf5_file):
    df = pd.read_csv(final_csv, low_memory=False)
    for _, target_dict in datasets.items():
        hdf5_file = f'{target_dict["main_path"]}/{target_dict["hdf5"]}'
        prefix = target_dict['prefix']
        df = replace_smiles_inchi_single(prefix, hdf5_file, df)

    # df = replace_smiles_inchi_total(final_hdf5_file, df)
    df = df.loc[:, ~df.columns.str.contains('Unnamed')]

    df.index+=1
    df.to_csv(final_csv)
    return df

def merge_hdf5(datasets, final_hdf5_file):
    with h5py.File(final_hdf5_file, 'w') as merged_file:
        for _, dataset_dict in datasets.items():
            hdf5_path = f'{dataset_dict["main_path"]}/{dataset_dict["hdf5"]}'
            with h5py.File(hdf5_path, 'r') as source_file:
                for data_i in source_file:
                    #if isinstance(source_file[data_i], h5py.Dataset):
                    new_name = dataset_dict['prefix'] + data_i
                    source_file.copy(source_file[data_i], merged_file, name=new_name)



def main(datasets):
    final_csv = f"{check_path}/final_all.csv"
    final_hdf5_file = f"{check_path}/final_all.hdf5"

    choose = int(input('1 - generate_final_csv_each \n2 - merge_hdf5 \n3 - merge_final_csv \n4 - replace smiles and inchi in final_csv \n'))
    if choose == 1:
        for target in datasets.items():
            generate_final_csv_each(target)
    elif choose ==2:
        merge_hdf5(datasets, final_hdf5_file)
    elif choose ==3:
        merge_final_csv(datasets, final_csv)
    elif choose == 4:
        replace_smiles_inchi(datasets, final_csv, final_hdf5_file)


if __name__ == '__main__':
    global check_path
    check_path = '/data/home/zhuyf/dataset_work/database/checkOptedGeoms'

    datasets = {
        'pubchemqc9' : {
            'prefix': 'Ba',

            'main_path': '/share/home/limg/analysis_lmg',
            'gs_dir': '/share/home/limg/pubchemqc/Gaussian_1_10_CNOF',
            'es_dir': '/data/home/zhuyf/pubchem_1_10_CNOF/Gaussian_1_10_CNOF_gau_excited_wb97xd',

            'csv_file': '9_CNOF_g_s_2.csv',
            'gs_error': 'pub_10_CNOF_error_list.csv',
            'es_error':'error_list_.csv',

            'hdf5':'9_CNOF_pub_molecules.hdf5',
        },
        'pubchemqc10' : {
            'prefix': 'Bb',
            'main_path': '/share/home/limg/analysis_lmg',
            'gs_dir': '/share/home/limg/pubchemqc/Gaussian_1_10_CNOF',
            'es_dir': '/data/home/zhuyf/pubchem_1_10_CNOF/Gaussian_1_10_CNOF_gau_excited_wb97xd',

            'csv_file': '10_CNOF_g_s_2.csv',
            'gs_error':'pub_10_CNOF_error_list.csv',
            'es_error':'error_list_.csv',

            'hdf5':'10_CNOF_pub_molecules.hdf5',
        },
        'qm9' : {
            'prefix': 'Aa',

            'main_path': '/data/home/zhuyf/dataset_work/database/qm9_database',
            'gs_dir': '/data/home/zhuyf/dataset_work/database/qm9_database/Gaussian',
            'es_dir': '/data/home/zhuyf/dataset_work/database/qm9_database/Gaussian_gau_excited_wb97xd',

            'csv_file': 'double_and_triple.csv',
            'gs_error':'error_list_double_and_triple.csv',
            'es_error': 'error_list_double_and_triple.csv',

            'hdf5':'qm9_props.hdf5',
        },
        'gdb11' : {
            'prefix': 'Ab',

            'main_path': '/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation',
            'gs_dir': '/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation/calculation_/gdb11_Gaussian_',
            'es_dir': '/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation/calculation_excited/gdb11_Gaussian_excited',

            'csv_file': 'new_mbkmeans_final_sample_1w.csv',
            'gs_error':'error_list_.csv',
            'es_error': 'error_list_excited.csv',

            'hdf5':'gdb11_props.hdf5',
        },
    }

    np.save(f'{check_path}/datasets_dict.npy', datasets)
    main(datasets)




