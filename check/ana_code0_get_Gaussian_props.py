import os
import numpy as np
import pandas as pd
import h5py
import json
import tempfile
import subprocess as sub

from openbabel import pybel
from rdkit import Chem

from sub_code_get_GS_opt_props import main_xyz as get_GS_opt_xyz

def write_molecule_data_to_hdf5(f, ind, smiles1, inchi1, smiles2,inchi2):
    molecule_group = f[ind]
    molecule_group.attrs['SMILES_PYBEL'] = smiles1
    molecule_group.attrs['INCHI_PYBEL'] = inchi1
    molecule_group.attrs['SMILES_RDKIT'] = smiles2
    molecule_group.attrs['INCHI_RDKIT'] = inchi2

def xyz_to_SMILESandINCHI_pybel(input_xyz_file):
    molecule = pybel.readfile("xyz", input_xyz_file).__next__()
    smi = molecule.write("smiles")
    inchi = molecule.write("inchi")
    return smi.split()[0], inchi.split()[0]

def xyz_to_SMILESandINCHI_rdkit(xyz_file):
    with tempfile.NamedTemporaryFile(delete=False) as temp_file_mol:
        mol_file = temp_file_mol.name
    sub.run(f'obabel -ixyz {xyz_file} -omol -O {mol_file}', shell=True)
    mol = Chem.MolFromMolFile(mol_file)
    smiles = Chem.MolToSmiles(mol)
    inchi = Chem.MolToInchi(mol)
    os.remove(mol_file)
    return smiles.strip(), inchi.strip()


def get_props_from_log(input_df,hdf5_file):
    err_list = []
    with h5py.File(hdf5_file, 'a') as f:
        for ind in input_df:
            print(ind)
            gau_dict = {
            'Aa':'/data/home/zhuyf/dataset_work/database/qm9_database/Gaussian',
            'Ab':'/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation/calculation_/gdb11_Gaussian_',
            'Ba':'/share/home/limg/pubchemqc/Gaussian_1_10_CNOF',
            'Bb':'/share/home/limg/pubchemqc/Gaussian_1_10_CNOF',
            }
            gau_path_GS = gau_dict[ind[:2]]
            log_file = f'{gau_path_GS}/{ind[2:]}/{ind[2:]}.log'
            with tempfile.NamedTemporaryFile(delete=False) as temp_file:
                xyz_file = temp_file.name
            get_GS_opt_xyz(log_file, xyz_file, work_type='opt+freq')
            try:
                smi1, inchi1 = xyz_to_SMILESandINCHI_pybel(xyz_file)
                smi2, inchi2 = xyz_to_SMILESandINCHI_rdkit(xyz_file)
                os.remove(xyz_file)
                write_molecule_data_to_hdf5(f, ind, smi1, inchi1, smi2, inchi2)
            except:
                err_list.append(ind)

    print(f'error_{err_list}')



if __name__ == '__main__':
    global label, gau_path_GS,  gau_path_ES, work_type

    work_type = 'opt+freq'


    main_path = f'/data/home/zhuyf/dataset_work/database/checkOptedGeoms'
    hdf5_file = f'{main_path}/final_all.hdf5'

    input_file = f'{main_path}/final_all.csv'

    input_df = pd.read_csv(input_file)['Index']


    get_props_from_log(input_df, hdf5_file)
    print('Done!')
