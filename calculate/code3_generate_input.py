'''
Author: zyf
'''

import numpy as np 
import os 
from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem import Draw
import pandas as pd
import shutil
from rdkit.Chem import Descriptors
# from openbabel import openbabel as ob


def generate_input(input_file):
    df=pd.read_csv(input_file)
    nums=df['Index']
    smiles=df['Smiles_GDB']
    for pair in zip(nums,smiles):
        
        mol = Chem.MolFromSmiles(pair[1])
        # 计算分子的自旋多重度
        num_radicals = Descriptors.NumRadicalElectrons(mol)
        spin_multiplicity = 2 * num_radicals + 1
        # 计算分子的总电荷
        total_charge = Chem.GetFormalCharge(mol)
        
        xyz_file = eval(choose['xyz_file'].replace('NUM', str(pair[0])))
        new_gjf_path = f'{gau_path}/{pair[0]}'
        os.makedirs(new_gjf_path,exist_ok=True)
        with open(xyz_file, 'r') as f1, open(f'{new_gjf_path}/{pair[0]}.{choose["suffix"]}', 'w') as f2:
            setting = choose['setting']
            f2.write('\n'.join(setting).replace('NUM', str(pair[0])))
            f2.write(f'{choose["cal_setting"]}{total_charge}   {spin_multiplicity}')
            if choose["xyz_need"] == True:
                lines=f1.readlines()
                num_atom = int(lines[0].strip())
                lines_coord = lines[2:2+num_atom]
                coord=[line.strip().split()[0:4] for line in lines_coord]
                xyz_string='\n'
                for atom in coord:
                    atom_symbol = atom[0]
                    atom_coords = '\t'.join(atom[1:]).replace('*^', 'e')
                    xyz_string += f"{atom_symbol}\t{atom_coords}\n"
                f2.write(xyz_string)
            end_line = '\n'.join(choose['end_setting'])
            f2.write(end_line)
            
        '''  
        if choose['qsub_file']:
            qsub_mode(pair[0])
        '''  

# def generate_input(input_file):
#     os.makedirs(f'./{label}_{choose_name}_{method}', exist_ok=True)
#     df=pd.read_csv(input_file)
#     nums=df['Index']
#     smiles=df['Smiles_GDB']
#     for i in nums:
#         input_xyz_file = f'./database/xyz_file_all/{i}.xyz' 
#         output_file = f'./{label}_{choose_name}_{method}/{pair[0]}.{choose["suffix"]}
#         xyz_to_gjf(input_xyz_file, output_mol_file)

'''  
def qsub_mode(num):
    os.makedirs(f'./{label}_{choose_name}_{method}/{num}', exist_ok=True)
    shutil.move(f'./{label}_{choose_name}_{method}/{num}.{choose["suffix"]}', f'./{choose_name}_{method}/{num}/{num}.{choose["suffix"]}')
    with open(choose['qsub_file'], 'r') as f1, open(f'./{label}_{choose_name}_{method}/qsub_{num}.sh', 'w') as f2:
        for line in f1:
            if 'input_name' in line:
                f2.write(f'{line.replace("input_name", str(num))}')
            else:
                f2.write(line)
'''

# def xyz_to_gjf(input_xyz_file, output_mol_file):
#     obConversion = ob.OBConversion()
#     obMolecule = ob.OBMol()

#     obConversion.SetInAndOutFormats("xyz", "gjf")

#     if obConversion.ReadFile(obMolecule, input_xyz_file):
#          obConversion.WriteFile(obMolecule, output_mol_file)
    
     

if __name__ == '__main__':
    #label = 'aromatic'
    label = 'bond_but_nonaromatic'
    
    global main_path,gau_path
    main_path = f'/data/home/zhuyf/dataset_work/database/qm9_database'
    old_gua_path = f'{main_path}/Gaussian'
    
    input_file=f'{main_path}/{label}.csv'
    
    global method, choose_name, choose
        
    state = input('\n1 - Ground-state (default: b3lyp+gd3bj/6-31G(d))\n2 - Excited-state (default: wb97xd/6-31G(d))\n\n')
    
    if state == '1':
        method='b3lyp'
        choose_name = 'gau'


    elif state == '2':
        method='wb97xd'
        choose_name = 'gau_excited'
    choose_dict = {'gau':{'setting':['%chk=NUM.chk', '%mem=8GB','%nprocshared=2', \
                    f'#P opt freq {method}/6-31G* empiricaldispersion=gd3bj SCF=VeryTight\n',\
                    'Title Card Required\n\n'], \
                    'cal_setting': '', \
                    'suffix':'gjf', \
                    'end_setting':['\n'], \
                    'qsub_file' : [],\
                    'xyz_file':"f'{main_path}/database/dsgdb9nsd_{NUM:06d}.xyz'",\
                    'xyz_need' : True
                    },
                   
                    'gau_excited_read':{'setting':[f'%oldchk={old_gua_path}/NUM/NUM.chk','%chk=NUM.chk', '%nprocshared=2', \
                    f'#P td=(nstates=10, 50-50)  {method}/6-31g(d) guess=read geom=check\n',\
                    'Title Card Required\n\n'], \
                    'cal_setting': '', \
                    'suffix':'gjf', \
                    'end_setting':['\n',''], \
                    'qsub_file' : [],
                    'xyz_file': "f'{main_path}/Gaussian/xyz_file_all/{NUM}.xyz'",\
                    'xyz_need' : False,
                    },
                    
                    'gau_excited':{'setting':['%chk=NUM.chk', '%mem=8GB', '%nprocshared=2', \
                    f'#P td=(nstates=10, 50-50)  {method}/6-31g(d)\n',\
                    'Title Card Required\n\n'], \
                    'cal_setting': '', \
                    'suffix':'gjf', \
                    'end_setting':['\n',''], \
                    'qsub_file' : [],
                    'xyz_file': "f'{main_path}/Gaussian/xyz_file_all/{NUM}.xyz'",\
                    'xyz_need' : True,
                    },
            }
    
    gau_path = f'{main_path}/Gaussian_{choose_name}_{method}'
    os.makedirs(gau_path, exist_ok=True)

    choose = choose_dict[choose_name]
    generate_input(input_file)

        
    print('DONE!')

