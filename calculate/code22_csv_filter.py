import pandas as pd
from code2_check_smiles import mv_image_and_tar 

#1
set1_file = './double_and_triple.csv' 
set2_file = './aromatic.csv'

#2
# set1_file = './smiles.csv' 
# set2_file = './double_and_triple.csv'
# set3_file = './aromatic.csv'

set1 = pd.read_csv(set1_file)
set2 = pd.read_csv(set2_file)

#2
# set3 = pd.read_csv(set3_file)

df3 = pd.DataFrame(columns=['Index', 'Smiles_GDB', 'Smiles_relaxed','InChI_GDB', 'InChI_relaxed', 'Image'])
#1
list_target = list(set(set1['Index']) - set(set2['Index']))

#2
# union_result = list(set(set2['Index']).union(set(set3['Index'])))
# list_target = list(set(set1['Index']) - set(union_result))

for i in list_target:
    df3=df3.append(set1[set1['Index'] == i], ignore_index=True)

df3.index+=1
df3 = df3.drop(df3.columns[-1], axis=1)

#1
df3.to_csv(f'./bond_but_nonaromatic.csv')
# mv_image_and_tar('bond_but_nonaromatic')

#2
# df3.to_csv(f'./single_bond_and_nonaromatic.csv')
# mv_image_and_tar('single_bond_and_nonaromatic')
