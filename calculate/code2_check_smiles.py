from rdkit import Chem 
import pandas as pd
import os 
import shutil
import tarfile
def check_smiles(choose):

    df=pd.read_csv(source_file)
    smis1=df['Smiles_GDB']
    index_list=list(range(1,len(smis1)+1))
    non_match_list=[]
    
    choose_dict = {1:'aromatic', 2:'double_and_triple'}
    
    df1 = pd.DataFrame(columns=df.columns)
    df2 = pd.DataFrame(columns=df.columns)

    for pair in zip(index_list, smis1):
        print(pair[0])
        mol = Chem.MolFromSmiles(pair[1])
        if (choose == 1 or choose ==3) and check_aromatic(mol):
            #true = check_aromatic(mol)
            #if true == True:
            df1=df1.append(df[df['Index'] == pair[0]], ignore_index=True)

        elif if (choose == 2 or choose == 3) and check_double_and_triple_bond(mol):
            #true = check_double_and_triple_bond(mol)
            #if true == True:
            df2=df2.append(df[df['Index'] == pair[0]], ignore_index=True)
     if choose == 1 or choose == 2:
        target_file_name = prefix+choose_dict[int(choose)]
        df_to_save = df1 if choose == 1 else df2
        df_to_save.index += 1
        df_to_save = df_to_save.drop(df_to_save.columns[-1], axis=1)
        # df_to_save = df_to_save.drop([col for col in df.columns if 'Unnamed' in col], axis=1)
        df_to_save.to_csv(f'{main_path}/{target_file_name}.csv')   
    elif choose ==3:
        target_file_name1 = prefix+choose_dict[1]
        target_file_name2 = prefix+choose_dict[2]
        df1.index+=1
        df1 = df1.drop(df1.columns[-1], axis=1)
        df1.to_csv(f'{main_path}/{target_file_name1}.csv')
        df2.index+=1
        df2 = df2.drop(df2.columns[-1], axis=1)
        df2.to_csv(f'{main_path}/{target_file_namecsv')

'''
def check_smiles(choose):

    df=pd.read_csv(source_file)
    smis1=df['Smiles_GDB']
    # smis2=df['smiles_relaxed']
    # inchi1=df['InChI_GDB']
    index_list=list(range(1,len(smis1)+1))
    # index_list=list(range(1,len(inchi1)+1))
    non_match_list=[]

    df2 = pd.DataFrame(columns=['Index', 'Smiles_GDB', 'Smiles_relaxed','InChI_GDB', 'InChI_relaxed'])

    for pair in zip(index_list, smis1):
    # for pair in zip(index_list, inchi1):
        print(pair[0])
        mol = Chem.MolFromSmiles(pair[1])
        # mol = Chem.MolFromInchi(pair[1])
        if choose == 1:
            true = check_aromatic(mol)

        elif choose == 2:
            true = check_double_and_triple_bond(mol)

        if true == True:
            # df2=df2.append(df[df['Smiles_GDB'] == pair[1]], ignore_index=True)
            df2=df2.append(df[df['Index'] == pair[0]], ignore_index=True)
            
    df2.index+=1
    df2 = df2.drop(df2.columns[-1], axis=1)
    df2.to_csv(f'./{target_file_name}.csv')
'''

def check_aromatic(mol):
    smarts = '[a]' 
    pattern = Chem.MolFromSmarts(smarts)
    is_aromatic = mol.HasSubstructMatch(pattern)

    if is_aromatic:
        # df2=df2.append(df.iloc[num], ignore_index=True)
        return True

def check_double_and_triple_bond(mol):

    # has_double_bond = any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in mol.GetBonds())
    # has_triple_bond = any(bond.GetBondType() == Chem.BondType.TRIPLE for bond in mol.GetBonds())
    if Chem.Mol.GetNumBonds(mol, Chem.BondType.DOUBLE) > 0 or Chem.Mol.GetNumBonds(mol, Chem.BondType.TRIPLE) > 0:
        # df2=df2.append(df.iloc[num], ignore_index=True)
        return True


def mv_image_and_tar(target_file_name):
    
    source_path = f'./{target_file_name}_images'
    os.makedirs(source_path, exist_ok=True)
    
    df2 = pd.read_csv(f'./{target_file_name}.csv')
    for i in df2.iloc[:, 1]:
        shutil.copyfile(f'./images/temp_{i}.png', f'./{target_file_name}_images/temp_{i}.png')
        
    target_file = f'./{target_file_name}_images.tar.gz' 
    with tarfile.open(target_file, 'w:gz') as tar:
        tar.add(source_path, arcname='')
    



    
if __name__ == '__main__':
    global source_file, target_file_name
    choose = input('1 -- aromatic\n2 -- double and triple bond\n3 - all\n')
    source_file = './smiles.csv'
    
    check_smiles(int(choose))
    # mv_image_and_tar(target_file_name)
