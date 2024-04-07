import numpy as np
import os
import pandas as pd
from code51_inchi_based_after_opt_check import xyz_to_inchi

def generate_xyz():

    for i in df:
        print(i)
        with open(f'./database/dsgdb9nsd_{int(i):06d}.xyz', 'r') as f1, open(f'{xyz_path}/{i}.xyz', 'w') as f2:
            
            lines=f1.readlines()[2:-3]
            coord=[line.strip().split()[0:4] for line in lines]
            xyz_string=''
            for atom in coord:
                atom_symbol = atom[0]
                atom_coords = '\t'.join(atom[1:]).replace('*^', 'e')
                xyz_string += f"{atom_symbol}\t{atom_coords}\n"
            num_atom = len(coord)
            f2.write(f'{num_atom}\ndsgdb9nsd_{int(i):06d}.xyz\n')
            f2.write(xyz_string)
            f2.write('\n')

def get_inchi(df):
    df2 = pd.DataFrame(columns=['Index', 'InChI_new'])
    for index_i in df:
        print(index_i)
        xyz_file = f'{xyz_path}/{index_i}.xyz'
        inchi = xyz_to_inchi(xyz_file)
        df2=df2.append({"Index": index_i, "InChI": inchi}, ignore_index=True)
    df2.index+=1
    df2.to_csv(f'{main_path}/old_inchi.csv')

if __name__ == '__main__':
    global main_path, xyz_path
    main_path = '/data/home/zhuyf/dataset_work/database/qm9_database'
    xyz_path=f'{main_path}/database/xyz_file_all'
    os.makedirs(xyz_path, exist_ok = True)
    input_file = f'{main_path}/smiles.csv'
    
    df=pd.read_csv(input_file)['Index']

    #generate_xyz(df)
    get_inchi(df)
    
