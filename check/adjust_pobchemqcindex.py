import pandas as pd

# 步骤1: 读取CSV文件
df = pd.read_csv('final_all.csv',low_memory=False)

# 步骤2: 找到Index列以"Ba"开头的行，并进行转换
# 对于满足条件的行，使用lambda函数将数字部分填充为9位
df['Index'] = df['Index'].apply(lambda x: f"Ba{x[2:].zfill(9)}" if x.startswith('Ba') else x)
df['Index'] = df['Index'].apply(lambda x: f"Bb{x[2:].zfill(9)}" if x.startswith('Bb') else x)

# 步骤4: 将修改后的DataFrame保存回CSV文件
df.to_csv('your_modified_file.csv', index=False)

