import os
import numpy as np
import networkx as nx
import networkx.algorithms.isomorphism as iso
from networkx.algorithms.isomorphism import categorical_node_match
import time
from ase.data import chemical_symbols
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def get_positions(filename):
    labels = []
    coords = []
    with open(filename, 'r') as file:
        file.readline()
        file.readline()
        for line in file:
            try:
                label, x, y, z =line.split()
            except:
                break
            labels.append(label)
            coords.append([x,y,z])
    return labels, coords



def set_graph(symbols, positions, MARGIN):
    # Angstrom unit
    Standard_bond_distances = {'H': {'H': 0.74, 'C':1.09, 'N':1.01, 'O':0.96, 'F':0.92} , \
                                'C': {'H': 1.09, 'C':1.54, 'N':1.47, 'O':1.43, 'F':1.35} , \
                                'N': {'H': 1.01, 'C':1.47, 'N':1.45, 'O':1.40, 'F':1.36} , \
                                'O': {'H': 0.96, 'C':1.43, 'N':1.40, 'O':1.48, 'F':1.42} , \
                                'F': {'H': 0.92, 'C':1.35, 'N':1.36, 'O':1.42, 'F':0.71} , \
                                }

    assert(len(symbols)==len(positions))
    size = len(symbols)
    cal_distance = lambda p1, p2: np.sqrt( (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2 )
    G=nx.Graph()
    for i, symbol in enumerate(symbols):
        G.add_node(i, symbol = symbol)

    for i in range(size):
        for j in range(i+1,size):
            # print(f'{cal_distance(positions[i],positions[j])} ---- {Standard_bond_distances[symbols[i]][symbols[j]]* MARGIN}')
            if(cal_distance(positions[i],positions[j]) < Standard_bond_distances[symbols[i]][symbols[j]]* MARGIN):
                G.add_edge(i,j)
    return G

def main_import(input_list, main_path, dataset_info):
    MARGIN = 1.15

    list_uniq_graph = []
    list_uniq_index = []
    list_duplicates = []
    graph_path=f'{main_path}/opted_graphs'
    os.makedirs(graph_path, exist_ok =True)
    prefix_gaupath={}
    for _, info_dict in dataset_info.items():
        prefix_gaupath[info_dict['prefix']] = info_dict['gs_dir']
    for ind_pre in input_list:

        prefix = ind_pre[:2]
        ind = ind_pre[2:]
        gau_path = prefix_gaupath[prefix]

        xyz_file=f'{gau_path}/{ind}/{ind}.xyz'
        symbols, positions = get_positions(xyz_file)
        positions = [[float(item) for item in sublist] for sublist in positions]
        G = set_graph(symbols, positions, MARGIN)
        save_graph(ind_pre, G, graph_path,symbols, positions)

        for _, graph in zip(list_uniq_index[::-1], list_uniq_graph[::-1]):
            if(nx.is_isomorphic(graph, G, node_match = categorical_node_match('symbol','X')  )):
                # print (_, ind)
                list_duplicates.append(ind_pre)
                break
        else:
            list_uniq_graph.append(G)
            list_uniq_index.append(ind)
    # print(list_duplicates)
    return list_duplicates
    # results = list(set(input_df) - set(list_duplicates))
    # return results


def main():
    global main_path, label, gau_path,input_file

    label = 'aromatic'
    # label = 'bond_but_nonaromatic'
    main_path = f'/data/home/zhuyf/dataset_work/database/qm9_database'
    gau_path = f'{main_path}/Gaussian'
    input_file = f'{main_path}/{label}.csv'
    error_list_file = f'{gau_path}/error_list_{label}.dat'
    error_list = pd.read_csv(error_list_file)

    input_df = pd.read_csv(input_file)
    input_df = input_df[~input_df['Index'].isin(error_list)]['Index']
    main_import(input_df, main_path)

def save_graph(ind, graph_i,graph_path,symbols, positions):
    pos = {i: positions[i][:2] for i in range(len(symbols))} # 使用原子的 x 和 y 坐标
    labels = {i: symbols[i] for i in range(len(symbols))}
    plt.figure()
    nx.draw(graph_i, pos, labels=labels, with_labels=True, node_size=800, node_color='lightblue')

    # # 添加边的标签
    # for (u, v, d) in graph_i.edges(data=True):
    #     x0, y0 = pos[u]
    #     x1, y1 = pos[v]
    #     label = f"{d['weight']:.2f}"  # 格式化标签到小数点后两位
    #     plt.text((x0 + x1) / 2, (y0 + y1) / 2, label, color='red')

    plt.savefig(f'{graph_path}/{ind}.png')
    plt.close()

def test():
    MARGIN = 1.15

    # 创建一个示例图
    symbols = ['C', 'C', 'O', 'H', 'H', 'H', 'H']
    positions = [
        [1.170635e+00, -1.486950e-01, -8.000000e-06],
        [-2.337860e-01, 3.996670e-01, -3.000000e-06],
        [-1.237575e+00, -2.768890e-01, 0.000000e+00],
        [1.154453e+00, -1.240792e+00, -3.230000e-04],
        [1.713039e+00, 2.194140e-01, -8.807590e-01],
        [1.712746e+00, 2.188760e-01, 8.811540e-01],
        [-3.007390e-01, 1.511780e+00, -9.000000e-06]
    ]

    # 你需要定义 Standard_bond_distances 和 MARGIN
    # Standard_bond_distances 是一个字典，包含每对原子符号对应的标准键长
    # MARGIN 是一个数值，用于调整判断两个原子是否连接的距离阈值
    # 例如：Standard_bond_distances = {'C': {'C': 1.54, 'O': 1.43, 'H': 1.09}, 'O': {'C': 1.43, 'O': 1.48, 'H': 0.96}, 'H': {'C': 1.09, 'O': 0.96, 'H': 0.74}}
    # MARGIN = 1.2

    G = set_graph(symbols, positions, MARGIN)

    # 绘制图
    pos = {i: positions[i][:2] for i in range(len(symbols))} # 使用原子的 x 和 y 坐标
    labels = {i: symbols[i] for i in range(len(symbols))}
    nx.draw(G, pos, labels=labels, with_labels=True, node_size=800, node_color='lightblue')

    # 显示图形
    plt.show()


if __name__ == '__main__':
    main()
