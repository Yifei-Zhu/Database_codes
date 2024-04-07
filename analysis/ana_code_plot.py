import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np

def count_rings():
    total_counts = pd.Series(dtype=str)
    for nu in range(2,11):
        #file_path = f'final_all_remove_duplicate_{nu}_atoms.csv'
        file_path = f'{label}_filtered_final_all_remove_duplicate_{nu}_atoms.csv'
        df = pd.read_csv(file_path)
        counts = df.groupby('RingNumber').size()
        total_counts = total_counts.add(counts, fill_value=0)
    print(total_counts)

    with open(f'{code_path}/ring_num.json', 'w') as json_file:
        json.dump(total_counts.to_dict(), json_file)
    #plot_ring_num(total_counts)
    plot_ring_num_ratio(total_counts)


def plot_ring_num_ratio(total_counts):
    # 对总数进行排序
    #total_counts.sort_values(ascending=True, inplace=True)
    total_sum = total_counts.sum()
    
    # 计算每个条目的百分比并绘制条形图
    bars = plt.bar(total_counts.index, 100 * total_counts / total_sum, color='#008E9B')
    plt.xlabel('Number of Rings', fontname='Arial', fontsize=18, fontweight='bold')
    plt.ylabel('Percentage of Molecules (%)', fontname='Arial', fontsize=18, fontweight='bold')

    # 为每个条形图添加百分比注释
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2., height,
                 f'{height:.3f}%', ha='center', va='bottom', fontsize=9)

    # 设置y轴的上限，使得有足够的空间显示注释
    plt.ylim(0, max(100 * total_counts / total_sum) * 1.1)
    plt.xticks(total_counts.index, rotation = 360)

    # 保存图表，确保code_path和label变量已正确传入
    plt.savefig(f'{code_path}/{label}ring_num_percentage.svg')


def plot_ring_num(total_counts):
    ax = total_counts.plot(kind='bar', color='#008E9B')
    plt.xlabel('Number of Ring', fontname='Arial', fontsize=18, fontweight='bold')
    plt.ylabel('Molecule Counts', fontname='Arial', fontsize=18,fontweight='bold')
    # plt.title('Category Counts Across Multiple CSV Files')
    for p in ax.patches:
        ax.annotate(str(p.get_height()), (p.get_x() + p.get_width() / 2., p.get_height()),
                ha='center', va='center', xytext=(0, 11), textcoords='offset points')
    plt.ylim(0, max(total_counts) * 1.1)
    plt.xticks(rotation = 360)

    plt.savefig(f'{code_path}/{label}ring_num.svg')

def count_molecular_categories():
    total_counts = pd.Series(dtype=str)
    for nu in range(2,11):
        file_path = f'final_all_remove_duplicate_{nu}_atoms.csv'
        #file_path = f'{label}_filtered_final_all_remove_duplicate_{nu}_atoms.csv'
        df = pd.read_csv(file_path)
        counts = df.groupby('CompoundType').size()
        total_counts = total_counts.add(counts, fill_value=0)
    with open(f'{code_path}/molecular_categories.json', 'w') as json_file:
        json.dump(total_counts.to_dict(), json_file)
    #plot_categories(total_counts)
    print(total_counts)
    plot_categories_ratio(total_counts)
    #plot_categories_pie(counts)

def plot_categories_ratio(total_counts):
    #total_counts.sort_values(inplace=True)
    
    # 计算总数以便计算每个条目的百分比
    total_sum = total_counts.sum()

    # plt.figure(figsize=(11, len(total_counts) / 2))  # 根据需要调整大小

    bars = plt.barh(total_counts.index, 100 * total_counts.values / total_sum, color='#008E9B')

    # 在条形图上方添加百分比数值
    for bar in bars:
        plt.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
                f'{bar.get_width():.1f}%', va='center')

    plt.xlabel('Percentage of Total Molecules', fontsize=18, fontweight='bold')
    plt.ylabel('Molecular Categories', fontsize=18, fontweight='bold')
    plt.xlim(0, max(100 * total_counts.values / total_sum) * 1.2)  # 调整为百分比的最大值的120%
    plt.title('Database Contents as Function of Molecular Categories')
    plt.tight_layout()
    plt.savefig(f'{code_path}/{label}molecular_categories_percentage.svg')


def plot_categories(total_counts):
    total_counts.sort_values(inplace=True)

    # plt.figure(figsize=(11, len(total_counts) / 2))  # Adjust the size as needed

    bars = plt.barh(total_counts.index, total_counts.values,color='#008E9B')

    # Add numbers on top of the bars
    for bar in bars:
        plt.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
                f'{int(bar.get_width())}', va='center')

    plt.xlabel('Molecule Counts', fontsize=18, fontweight='bold')
    plt.ylabel('Molecular Categories', fontsize=18, fontweight='bold')
    plt.xlim(0, max(total_counts) * 1.2)
    plt.title('Database contents as function of molecular categories')
    plt.tight_layout()
    plt.savefig(f'{code_path}/{label}molecular_categories.svg')

def plot_categories_pie(counts):
    # academic_colors = ['#008F7A','#0081CF','#845EC2', '#D65DB1', '#FF6F91', '#FF9671', '#FFC75F', '#F9F871']
    # plt.figure(figsize=(10, 10))
    academic_colors = ['#F90C0C','#FFEBCD', '#00A0BA','#00BDB4','#4AD69E','#A3EB83','#616EA5','#ED992A',]
    explode = [0.2 if value < 5 else 0 for value in counts.values / counts.sum() * 100]
    # explode = [0.1 if value == max(counts.values) else 0 for value in counts.values]
    # explode=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,]
    # explode=[0,0,0,0,0,0,0,0]
    plt.pie(counts.values, labels=counts.index, autopct='%1.1f%%', startangle=140, explode=explode,colors=academic_colors[:len(counts)])

    # # Set the figure size
    # plt.figure(figsize=(12, 12))  # Adjust the size as needed

    # # Plotting the pie chart
    # wedges, texts = plt.pie(counts.values, labels=counts.index, startangle=140)

    # # Increase the font size for the labels
    # for text in texts:
    #     text.set_size('large')

    # # Place labels with lines and increase font size
    # for i, p in enumerate(wedges):
    #     ang = (p.theta2 - p.theta1) / 2 + p.theta1
    #     y = np.sin(np.deg2rad(ang))
    #     x = np.cos(np.deg2rad(ang))
    #     connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    #     kw = dict(arrowprops=dict(arrowstyle="->", connectionstyle=connectionstyle), zorder=0, va="center")
    #     plt.annotate(f'{counts.values[i]}', xy=(x, y), xytext=(1.2*x, 1.2*y),
    #                 horizontalalignment='center', **kw)

    plt.title('Molecular Category Distribution')
    plt.savefig(f'{code_path}/{label}molecular_categories.svg')




if __name__ == '__main__':

    from matplotlib import rcParams
    
    
    config = {
        "font.family":'Times New Roman',  # 设置字体类型
        "axes.unicode_minus": False, #解决负号无法显示的问题
        "axes.labelsize":18,
    }
    rcParams.update(config)
                        
  
    global label
    label = ''


    code_path = '/data/home/zhuyf/dataset_work/database/analysis'
    main_path = '/data/home/zhuyf/dataset_work/database/qm9_database'

    #count_rings()
    count_molecular_categories()
    
