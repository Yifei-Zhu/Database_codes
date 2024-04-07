import os
import pandas as pd

def ric_main(input_df):
    error_list = []
    for index_i in input_df:
        gau_file = f'{gau_path}/{index_i}/{index_i}.log'
        print(gau_file)
        s = compare_ric_single(gau_file)
        if s == False:
            error_list.append(index_i)
        # break
    print(f'Number of molecules with different initial and optimized redundant internal coordinates:{len(error_list)}')
    with open(f'{gau_path}/{label}_error_redundant_ic.dat', 'w') as f1:
        for item in error_list:
            f1.write(str(item) + '\n')

def compare_ric_single(gau_file):
    initial_ric = get_initial_ric(gau_file)
    opted_ric =get_opted_ric(gau_file)
    if initial_ric != opted_ric:
        print(gau_file)
        s = False
    else:
        s = True
        with open(f'{os.path.dirname(gau_file)}/redundant_ic.dat', 'w') as f1:
            for item in opted_ric:
                f1.write(str(item) + '\n')
    return s

def get_initial_ric(gau_file):
    file_dir_name = os.path.splitext(gau_file)[0]
    file_extension = os.path.splitext(gau_file)[1]
    if os.path.exists(f'{file_dir_name}_1{file_extension}'):
        gau_file=f'{file_dir_name}_1{file_extension}'
    gau_keyword = 'Initial Parameters'
    line_nu = get_target_line_nu(gau_file, gau_keyword, count=1)
    ric = get_redundant_iC(gau_file, line_nu)
    return ric

def get_opted_ric(gau_file):
    gau_keyword = 'Optimized Parameters'
    line_nu = get_target_line_nu(gau_file, gau_keyword, count=-1)
    ric = get_redundant_iC(gau_file, line_nu)
    return ric

def get_target_line_nu(gau_file, gau_keyword, count=None):
    with open(gau_file, 'r') as f1:
        lines = f1.readlines()[::count]
        line_nu = -1
        for line in lines:
            line_nu+=1
            if gau_keyword in line:
                break
    return line_nu*count + int(count+0.5)


def get_redundant_iC(gau_file, line_nu):
    ric = []
    with open(gau_file, 'r') as f1:
        lines = f1.readlines()[line_nu+4:]
        for line in lines:
            if '!' not in line:
                break
            ric.append(line.split()[1:3])
    return ric


def main():
    global main_path, label, gau_path

    #label = 'aromatic'
    #label = 'bond_but_nonaromatic'
    '''
    label = 'double_and_triple'
    main_path = f'/data/home/zhuyf/dataset_work/database/qm9_database'
    gau_path = f'{main_path}/Gaussian'
    '''

    label=''
    #label='excited'
    prefix = 'gdb11_'
    main_path = '/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation'
    gau_path = f'{main_path}/calculation_{label}/{prefix}Gaussian_{label}'
  
    input_file = f'{gau_path}/{label}_error_inchi.csv'
    error_list_file = f'{gau_path}/error_list_{label}.csv'
    error_list = pd.read_csv(error_list_file)['Index']
    input_df = pd.read_csv(input_file)
    input_df = input_df[~input_df['Index'].isin(error_list)]['Index']


    ric_main(input_df)

if __name__ == '__main__':
    main()
