import pandas as pd

# 读取CSV文件
df = pd.read_csv('final_all_remove_duplicate.csv')

# 如果'CompoundType'列存在，则删除它
if 'Unnamed: 0.1.1' in df.columns:
    df = df.drop(columns=['Unnamed: 0.1.1'])
if 'Unnamed: 0.1' in df.columns:
    df = df.drop(columns=['Unnamed: 0.1'])
if 'Unnamed: 0' in df.columns:
    df = df.drop(columns=['Unnamed: 0'])
df = df.drop(columns=['CompoundType'])


# 将'CoumpoundType'列重命名为'CompoundType'
df = df.rename(columns={'CoumpoundType': 'CompoundType'})

# 将修改后的DataFrame保存到新的CSV文件
df.to_csv('final_all_remove_duplicate.csv', index=False)



