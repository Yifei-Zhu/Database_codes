import pickle
import pandas as pd
import subprocess as sub
from collections import Counter

from FunctionalGroups.searchFunctionalGroups import main_searchGroupsImport as serch_FG
from rdkit.Chem import FragmentCatalog

from rdkit import Chem
from rdkit.Chem import Draw

def analyzeFGs_main(input_file):
    counter = Counter()
    none_list = []
    df = pd.read_csv(input_file, index_col=0)[col_label]
    error_list= []
    for smi in df:
        try:
            fgs_id = serch_FG(smi)
        except:
            error_list.append(smi)
            print(error_list)
            continue
        if fgs_id == None:
            none_list.append(smi)
        counter.update(fgs_id)
    all_dict = {
    "count_id":counter,
    "list": none_list
    }
    with open(f'{code_path}/FunctionalGroups_count.pkl', 'wb') as file:
        pickle.dump(all_dict, file)


def draw_img_of_top():
    with open(f'{code_path}/FunctionalGroups_count.pkl', 'rb') as file:
        df = pickle.load(file)
    fg_file=f'{code_path}/FunctionalGroups/FunctionalGroups_for_plot.txt'
    fparams = FragmentCatalog.FragCatParams(1,6,fg_file)
    mols=[]
    top=24
    sorted_dict = sorted(df['count_id'].items(), key=lambda x: x[1], reverse=True)[:top]
    for key, value in sorted_dict:
        funcgroup = fparams.GetFuncGroup(key)
        name=f"{funcgroup.GetProp('_Name')}_{key}"
        print(f"{name}: {value}")
        mols.append(funcgroup)
    img=Draw.MolsToGridImage(mols,molsPerRow=6)
    img.save(f'{code_path}/fgs_of_top_{top}.png')
    print(sorted(df['count_id'].keys()))
    print(len(df['count_id']))
    print(df['list'])

    # for i in [3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 17, 19, 24, 25, 26, 32, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 52, 53, 55, 57, 58, 59, 60, 62, 63, 64, 65, 67, 69, 71, 76, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 92, 93, 94, 97, 108, 111, 115, 119]:
    #     funcgroup = fparams.GetFuncGroup(i)
    #     name=f"{funcgroup.GetProp('_Name')}_{i}"
    #     print(name)


def checkmol_main(input_file):
    df = pd.read_csv(input_file, index_col=0)[col_label]
    # smiles_to_mol(df)
    run_checkmol(df)

def smiles_to_mol(df):
    for nu,smi in enumerate(df):
        mol = Chem.MolFromSmiles(smi)
        mol_block = Chem.MolToMolBlock(mol)
        mol_file_path = f'{code_path}/FunctionalGroups/checkmol/mols/index_{nu}.mol'
        with open(mol_file_path, 'w') as f:
            f.write(mol_block)

def run_checkmol(df):
    counter = Counter()
    for nu,_ in enumerate(df):
        result = sub.run(f'{code_path}/FunctionalGroups/checkmol/checkmol-0.5b-linux-x86_64 -e {code_path}/FunctionalGroups/checkmol/mols/index_{nu}.mol',shell=True, capture_output=True, text=True)
        out = result.stdout.strip().split('\n')
        counter.update(out)
    with open(f'{code_path}/FunctionalGroups_count_checkmol.pkl', 'wb') as file:
        pickle.dump(counter, file)

def checkmol_top():
    with open(f'{code_path}/FunctionalGroups_count_checkmol.pkl', 'rb') as file:
        df = pickle.load(file)
    #print(df)
    print(len(df))

if __name__ == '__main__':
    code_path = '/data/home/zhuyf/dataset_work/database/analysis'
    main_path = '/data/home/zhuyf/dataset_work/database/analysis'

    smi_file = f'{main_path}/final_all_remove_duplicate.csv'
    #smi_file = f'{main_path}/filtered_final_all_remove_duplicate.csv'
    #smi_file = f'{main_path}/B_filtered_final_all_remove_duplicate.csv'

    global col_label
    col_label='Smiles'

    analyzeFGs_main(smi_file)
    draw_img_of_top()

    # checkmol_main(smi_file)
    # checkmol_top()
