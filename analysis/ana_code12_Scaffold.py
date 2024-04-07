import os
import re
import numpy as np
import pandas as pd
import json
import pickle
from collections import defaultdict,Counter

import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold


def analyzeScaffold(smiles_list, inchi_list):
    mols=[]
    for smiles, inchi in zip(smiles_list, inchi_list):

        mol = None
        if smiles != '1':
            mol = Chem.MolFromSmiles(smiles)
        if mol is None:  # 如果SMILES转换失败
            mol = Chem.MolFromInchi(inchi)  # 尝试使用InChI字符串

        mols.append(mol)

    mol_scaffolds = [MurckoScaffold.GetScaffoldForMol(mol) for mol in mols if mol]
    # for i in drugbank_atomic_scaffolds:
    #     i.Compute2DCoords()
    if label == 'C':
    ### Atoms are treated as C Bonds are treated as single ones
        grafh_scaffolds = [MurckoScaffold.MakeScaffoldGeneric(s) for s in mol_scaffolds]
    ### Normal
    else:
        grafh_scaffolds = [MurckoScaffold.GetScaffoldForMol(s) for s in mol_scaffolds]
    scaffold_smiles = [Chem.MolToSmiles(scaffold) for scaffold in grafh_scaffolds if scaffold != None]

    lists=[]
    for nu,(i,j) in enumerate(zip(mols, scaffold_smiles)):
        if i. GetRingInfo().NumRings() == 0:
            if j != '':
                lists.append(smiles_list[nu])

    counter=Counter(scaffold_smiles)
    with open(f'{code_path}/{label2}MurckoScaffold_all_{label}.pkl', 'wb') as file:
        pickle.dump(counter, file)

    sorted_dict = dict(sorted(counter.items(), key=lambda item: item[1], reverse=True)[:20])
    for key, value in sorted_dict.items():
        print(f'{key}: {value}')
    len_functional_groups_count = len(counter)
    print(f"Length of the loaded object: {len_functional_groups_count}")


def analyzeScaffold_main(input_file):
    smiles_list = pd.read_csv(input_file, index_col=0)[col_label1]
    inchi_list = pd.read_csv(input_file, index_col=0)[col_label2]
    analyzeScaffold(smiles_list, inchi_list)

def get_scaffold_fig(df):

    sorted_dict = sorted(df.items(), key=lambda x: x[1], reverse=True)[:20]

    for key, value in sorted_dict:
        print(f"{key}: {value}")

    smis = [i[0] for i in sorted_dict]
    print(smis)
    molecules = [Chem.MolFromSmiles(smiles) for smiles in smis]

    # 生成分子图
    img = Draw.MolsToGridImage(molecules, molsPerRow=4, subImgSize=(200,200))
    img.save(f'{code_path}/scaffolds_img_{label}.png')


def pie_plot(data):
    # 筛选数值大于1000的数据并计算其他的总和
    filtered_data = {k: v for k, v in data.items() if v > 2000}
    filtered_data["Other"] = sum(v for k, v in data.items() if v <= 2000)

    # 数据准备
    labels = filtered_data.keys()
    sizes = filtered_data.values()

    # 绘制饼状图
    plt.figure(figsize=(10, 10))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, textprops={'fontsize': 16}, pctdistance=0.85)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.title('Pie Chart of Categories with Values > 2000 and Others')
    plt.tight_layout()
    plt.savefig(f'{code_path}/scaffold_pie_{label}.png')

def plot_main():
    with open(f'{code_path}/MurckoScaffold_all_{label}.pkl', 'rb') as file:
        df = pickle.load(file)
    pie_plot(df)
    get_scaffold_fig(df)

if __name__ == '__main__':
    code_path = '/data/home/zhuyf/dataset_work/database/analysis'
    main_path = '/data/home/zhuyf/dataset_work/database/analysis'

    global label, label2
    choose = int(input('(1) Normal or (2) treat all atoms as C and all bonds as single type\n'))
    # choose =1

    label='atom' if choose == 1 else 'C'
    
    label2=''
    smi_file = f'{main_path}/final_all_remove_duplicate.csv'
    # smi_file = f'{main_path}/{label2}_filtered_final_all_remove_duplicate.csv'
    # smi_file = f'{main_path}/a.csv'
    global col_label1, col_label2
    col_label1 = 'Smiles_rdkit'
    col_label2 = 'InchI_pybel'

    analyzeScaffold_main(smi_file)
    #plot_main()
