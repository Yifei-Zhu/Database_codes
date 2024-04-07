import numpy as np 
import os 
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pandas as pd
# import sys
# sys.path.append('/data/home/zhuyf/CODES/rdkit_tools')
# from draw_mols_from_smiles import draw_from_smiles

def split_smiles():
    num_file=len(os.listdir(f'{path}/database'))

    # lines_list1, lines_list2 = [], []
    # img_list = []
    
    df = pd.DataFrame(columns=['Index', 'Smiles_GDB', 'Smiles_relaxed','InChI_GDB', 'InChI_relaxed'])
    # df2 = pd.DataFrame(columns=['Index', 'Smiles_GDB', 'Smiles_relaxed','InChI_GDB', 'InChI_relaxed'])
    df_zwitterions = pd.DataFrame(columns=['Index', 'Smiles_GDB', 'Smiles_relaxed','InChI_GDB', 'InChI_relaxed'])
    
    for i in range(num_file):
        with open(f'{path}/database/dsgdb9nsd_{i+1:06d}.xyz', 'r') as file:
            lines = file.readlines()
            line_i=lines[-2].strip()
            line_j=lines[-1].strip()
            # lines_list1.append(line_i.split()[0])
            # lines_list2.append(line_i.split()[1])
            smiles_gdb=line_i.split()[0]
            if '+' in smiles_gdb:
                df_zwitterions.append({"Index": i+1, "Smiles_GDB": line_i.split()[0], "Smiles_relaxed": line_i.split()[1], "InChI_GDB": line_j.split()[0], "InChI_relaxed": line_j.split()[1]}, ignore_index=True)
                continue
            
            try:
                mol = Chem.MolFromSmiles(line_i.split()[0])
                # mol = Chem.MolFromInchi(line_j.split()[0])
                drawer = Draw.rdMolDraw2D.MolDraw2DCairo(300, 300)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                image_path = f'{path}/images/temp_{i+1}.png'
                drawer.WriteDrawingText(image_path)
            except:
                print(i+1)
                # df2 = df2.append({"Index": i+1, "Smiles_GDB": line_i.split()[0], "Smiles_relaxed": line_i.split()[1], "InChI_GDB": line_j.split()[0], "InChI_relaxed": line_j.split()[1]}, ignore_index=True)
                continue

            df = df.append({"Index": i+1, "Smiles_GDB": line_i.split()[0], "Smiles_relaxed": line_i.split()[1], "InChI_GDB": line_j.split()[0], "InChI_relaxed": line_j.split()[1]},ignore_index=True)

    df.index = df.index + 1
    df_zwitterions.index = df_zwitterions.index + 1

    df.to_csv(f'{path}/smiles.csv')
    df_zwitterions.to_csv('./zwitterions.csv')
    # df2.to_csv('./smiles_fail.csv')
    # df.to_excel('./smiles.xlsx')
   

if __name__ == '__main__':
    global path
    path = '/data/home/zhuyf/dataset_work/database/qm9_database'
    split_smiles()
