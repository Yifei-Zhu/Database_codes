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
from sub_code_get_GS_opt_props import main_import as get_GS_props
from sub_code_get_ES_opt_props import main_import as get_ES_props

def write_molecule_data_to_hdf5(f, ind, smiles1, inchi1, smiles2,inchi2,  properties_ground, properties_excited):
    """
    将一个分子的数据写入HDF5文件

    :param hdf5_file: HDF5文件对象或路径。
    :param cid: 分子的CID。
    :param smiles: 分子的SMILES码。
    :param properties_ground: 分子基态的性质。
    :param properties_excited: 分子激发态的性质。
    """
    # 为当前分子创建一个组
    try:
        del f[str(ind)]
    except:
        pass
    molecule_group = f.create_group(str(ind))
    molecule_group.attrs['SMILES_PYBEL'] = smiles1
    molecule_group.attrs['INCHI_PYBEL'] = inchi1
    molecule_group.attrs['SMILES_RDKIT'] = smiles2
    molecule_group.attrs['INCHI_RDKIT'] = inchi2

    # 分别为基态和激发态性质创建数据集
    for state_name, properties in [('ground_state', properties_ground), ('excited_state', properties_excited)]:
        state_group = molecule_group.create_group(state_name)
        for prop_name, prop_value in properties.items():

            if isinstance(prop_value, dict):
                prop_json = json.dumps(prop_value)
                dt = h5py.special_dtype(vlen=str)
                # 使用可变长度字符串数据类型来创建数据集
                dset = state_group.create_dataset(prop_name, (1,), dtype=dt, compression='lzf')
                dset[0] = prop_json
            else:
                state_group.create_dataset(prop_name, data=prop_value, compression='lzf')

def xyz_to_SMILESandINCHI_pybel(input_xyz_file):
    molecule = pybel.readfile("xyz", input_xyz_file).__next__()
    try:
        smi = molecule.write("smiles").split()[0]
    except:
        smi = '1'
    try:
        inchi = molecule.write("inchi").split()[0]
    except:
        inchi = '1'
    return smi, inchi

def xyz_to_SMILESandINCHI_rdkit(xyz_file):
    with tempfile.NamedTemporaryFile(delete=False) as temp_file_mol:
        mol_file = temp_file_mol.name
    sub.run(f'obabel -ixyz {xyz_file} -omol -O {mol_file}', shell=True)
    try:
        mol = Chem.MolFromMolFile(mol_file)
        smiles = Chem.MolToSmiles(mol).strip()
        inchi = Chem.MolToInchi(mol).strip()
    except:
        smiles = '1'
        inchi = '1'

    os.remove(mol_file)
    return smiles, inchi


def get_props_from_log(input_df,hdf5_file):
    err_list = []
    with h5py.File(hdf5_file, 'a') as f:
        for ind in input_df:
            ind = str(ind).zfill(9)
            log_file = f'{gau_path_GS}/{ind}/{ind}.log'
            with tempfile.NamedTemporaryFile(delete=False) as temp_file:
                xyz_file = temp_file.name
            get_GS_opt_xyz(log_file, xyz_file, work_type='opt+freq')
            # try:
            smi1, inchi1 = xyz_to_SMILESandINCHI_pybel(xyz_file)
            smi2, inchi2 = xyz_to_SMILESandINCHI_rdkit(xyz_file)
            # except:
            #     err_list.append(ind)


            os.remove(xyz_file)
            log_file_GS=f'{gau_path_GS}/{ind}/{ind}.log'
            log_file_ES=f'{gau_path_ES}/{ind}/{ind}.log'
            properties_ground = get_GS_props(log_file_GS, work_type1)
            properties_excited = get_ES_props(log_file_ES, work_type2)
            write_molecule_data_to_hdf5(f, ind, smi1, inchi1, smi2, inchi2, properties_ground, properties_excited)
    print(err_list)

if __name__ == '__main__':
    # global label, gau_path, state, work_type
    global label, gau_path_GS,  gau_path_ES, work_type

    label = ''

    work_type1 = 'opt+freq'
    work_type2 = 'sp'


    main_path = f'/share/home/limg/analysis_lmg'
    hdf5_file = f'{main_path}/10_CNOF_pub_molecules.hdf5'

    gau_path_GS = f'/share/home/limg/pubchemqc/Gaussian_1_10_CNOF'
    gau_path_ES = f'/data/home/zhuyf/pubchem_1_10_CNOF/Gaussian_1_10_CNOF_gau_excited_wb97xd'

    try:
        error_list_path = f'{gau_path_GS}/error_list_{label}.csv'
        error_list= pd.read_csv(error_list_path)['Index']
    except FileNotFoundError:
        error_list= []

    input_file = f'/share/home/limg/analysis_lmg/10_CNOF_g_s_2.csv'

    input_df = pd.read_csv(input_file)

    input_df = input_df[~input_df['Index'].isin(error_list)]['Index']

    # input_df=['010796708', '055285513', '055286488', '072207987', '072651574', '074895648']
    input_df=['005369015', '005850345', '010290872', '073411456', '073754019', '059303776', '059347742']



    get_props_from_log(input_df, hdf5_file)
    print('Done!')
