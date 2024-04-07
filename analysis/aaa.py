import pandas as pd

# 读取CSV文件
df = pd.read_csv('final_all_remove_duplicate.csv')

# 筛选出'Index'列前两位为'Aa'或'Ab'的行
filtered_df = df[df['Index'].str.startswith(('Aa', 'Ab'))]

# 将筛选后的DataFrame保存到新的CSV文件
filtered_df.to_csv('A_filtered_final_all_remove_duplicate.csv', index=False)



