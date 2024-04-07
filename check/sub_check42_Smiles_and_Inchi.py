import h5py
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import defaultdict
import os
import shutil
import json

def check_by_smiles(hdf5_file):
    uni_list = []
    # dup_list = []
    dup_name_list = []
    with h5py.File(hdf5_file, 'r') as file:
        # for dataset_name in file:
        dataset_names = list(file)
        for dataset_name in reversed(dataset_names):

            dataset = file[dataset_name]
            smiles =  dataset.attrs['SMILES_PYBEL']
            if smiles not in uni_list:
                uni_list.append(smiles)
            else:
                # dup_list.append(smiles)
                dup_name_list.append(dataset_name)
    print(len(dup_name_list))
    return dup_name_list

def check_by_inchi(hdf5_file):
    uni_list = []
    # dup_list = []
    dup_name_list = []
    with h5py.File(hdf5_file, 'r') as file:
        # for dataset_name in file:
        dataset_names = list(file)
        for dataset_name in reversed(dataset_names):
            dataset = file[dataset_name]
            inchi = dataset.attrs['INCHI_PYBEL']
            if inchi not in uni_list:
                uni_list.append(inchi)
            else:
                # dup_list.append(inchi)
                dup_name_list.append(dataset_name)
    print(len(dup_name_list))
    return dup_name_list


def get_duplicate_file_smiles(hdf5_file, dataset_info, json_file):
    with open(json_file, 'r') as jsonfile:
        dup_list = json.load(jsonfile)
    dup_smiles = []
    with h5py.File(hdf5_file, 'r') as file:
        for dataset_name in dup_list:
            dataset = file[dataset_name]
            smiles = dataset.attrs['SMILES']
            dup_smiles.append(smiles)

    os.makedirs('/data/home/zhuyf/dataset_work/database/checkOptedGeoms/duplicates_smiles',exist_ok=True)
    prefix_gaupath = {}
    for _, info_dict in dataset_info.items():
        prefix_gaupath[info_dict['prefix']] = info_dict['gs_dir']
    with h5py.File(hdf5_file, 'r') as file:
        for dataset_name in file:
            dataset = file[dataset_name]
            smiles = dataset.attrs['SMILES']
            if smiles in dup_smiles:
                prefix = dataset_name[:2]
                ind = dataset_name[2:]
                log_file = f'{prefix_gaupath[prefix]}/{ind}/{ind}.log'
                gjf_file = f'{prefix_gaupath[prefix]}/{ind}/{ind}.gjf'
                new_path = f'/data/home/zhuyf/dataset_work/database/checkOptedGeoms/duplicates_smiles/{smiles}'
                os.makedirs(new_path, exist_ok=True)
                shutil.copyfile(log_file, os.path.join(new_path, f'{ind}.log'))
                shutil.copyfile(gjf_file, os.path.join(new_path, f'{ind}.gjf'))



def get_duplicate_file_inchi(hdf5_file, dataset_info, json_file):
    with open(json_file, 'r') as jsonfile:
        dup_list = json.load(jsonfile)

    dup_inchi = []
    with h5py.File(hdf5_file, 'r') as file:
        for dataset_name in dup_list:
            dataset = file[dataset_name]
            inchi = dataset.attrs['INCHI']
            dup_inchi.append(inchi)
    os.makedirs('/data/home/zhuyf/dataset_work/database/checkOptedGeoms/duplicates_inchi',exist_ok=True)
    prefix_gaupath = {}
    for _, info_dict in dataset_info.items():
        prefix_gaupath[info_dict['prefix']] = info_dict['es_dir']
    with h5py.File(hdf5_file, 'r') as file:
        for dataset_name in file:

            dataset = file[dataset_name]
            inchi = dataset.attrs['INCHI']

            if inchi in dup_inchi:

                prefix = dataset_name[:2]
                ind = dataset_name[2:]
                gjf_file = f'{prefix_gaupath[prefix]}/{ind}/{ind}.gjf'
                new_path = f'/data/home/zhuyf/dataset_work/database/checkOptedGeoms/duplicates_inchi/{inchi.replace("/", "_")}'
                os.makedirs(new_path, exist_ok=True)
                shutil.copyfile(gjf_file, os.path.join(new_path, f'{dataset_name}.gjf'))


if __name__ == '__main__':
    hdf5_file = 'final_all.hdf5'
    uni_list, dup_list=check_by_smiles(hdf5_file)
    print('dup_list:')
    print(dup_list)
