# 假设 "FunctionalGroups_count.pkl" 文件已经存在于当前工作目录中，我们将使用 Python 来读取这个文件
# 并打印其长度。首先，我们需要加载必要的库。
import pickle

with open('MurckoScaffold_all_atom.pkl', 'rb') as file:
    functional_groups_count = pickle.load(file)

sorted_dict = dict(sorted(functional_groups_count.items(), key=lambda item: item[1], reverse=True)[:20])

for key, value in sorted_dict.items():
    print(f'{key}: {value}')

len_functional_groups_count = len(functional_groups_count)
print(f"Length of the loaded object: {len_functional_groups_count}")
