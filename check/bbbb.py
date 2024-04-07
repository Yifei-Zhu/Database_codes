import pandas as pd

# 读取CSV文件
df = pd.read_csv('final_all_remove_duplicate_rev.csv')

# 删除不需要的列
#df.drop(columns=['Unnamed: 0.1'], inplace=True)
df.drop(columns=['index1'], inplace=True)

# 将修改后的DataFrame保存到新的CSV文件
df.to_csv('final_all_remove_duplicate_rev.csv', index=False)
