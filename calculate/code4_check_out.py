import numpy as np 
import os 
import pandas as pd
import shutil
import json
import subprocess as sub
import re
import time

def check_status_gau(num):
    error_check=False
    file_path = f'{gau_path}/{num}/{num}.log'
    with open(file_path, "r") as file:
        lines = file.readlines()
        last_line = lines[-1].strip()
        if "Normal termination" not in last_line:
            print(num)
            error_check=True
    return error_check

def check_main(df,error_list_file):
    # error_list = []
    nums=df['Index']
    df2 = pd.DataFrame(columns=['Index','error_type'])
    for num in nums:
        print(num)
        error_check = check_status_gau(num)
        if error_check == True:
            df2.loc[len(df2)]=[num, '']
        # error_list = check_status_gau(num, error_list)
    df2.index+=1
    df2.to_csv(error_list_file)  
    # np.savetxt(error_list_file, np.array(error_list), fmt='%d')
    # print(f'Total error:{len(error_list)}')


def error_type(df, error_list_file):
    # error_dict = {}
    # error_list = np.loadtxt(error_list_file, dtype=str).tolist()     
    # if not isinstance(error_list, list):
    #     error_list = [error_list]
        
    df2=pd.read_csv(error_list_file)

    for err_nu in df2['Index']:
        with open(f'{gau_path}/{err_nu}/{err_nu}.log', 'r') as f1:
            print(f'{gau_path}/{err_nu}/{err_nu}.log')
            lines=f1.readlines()[-10:-1]
            find=False
            for line in lines:
                if 'Error termination via Lnk1e' in line:
                    find = True
                    err_type = os.path.splitext(os.path.basename(line.split()[5]))[0]
                    # error_dict[int(err_nu)]=err_type
                    df2.loc[df2['Index'] == int(err_nu), 'error_type'] = err_type
                    break
            # if find == False:
                # error_dict[int(err_nu)]=''

    df2.to_csv(error_list_file, index=False)
    # with open(error_list_file, 'w') as file:
    #     for key, value in error_dict.items():
    #         file.write(f"{key}: {value}\n")


def read_dict_from_file(file_path):
    dictionary = {}
    with open(file_path, "r") as file:
        lines = file.readlines()
        for line in lines:
            key, value = line.strip().split(": ")
            dictionary[key] = value
    return dictionary

'''
def process_error_(value_to_find):
    def decorator(func):
        def wrapper(dictionary_path, *args, **kwargs):
            dictionary = read_dict_from_file(dictionary_path)
            
            keys_with_value = []
            for key, val in dictionary.items():
                if val == value_to_find:
                    keys_with_value.append(key)
                    

            for key in keys_with_value:
                out = open(f'{path}/{key}/{key}_1.log', 'r')
                lines1 = out.readlines()
                lines1.reverse()
                lines = []
                for line in lines1:
                    if 'Initial Parameters' in line:
                        break
                    else:
                        lines.append(line)
                lines.reverse()

                # os.rename(f'{path}/{key}/{key}.gjf', f'{path}/{key}/{key}_1.gjf')
                # os.rename(f'{path}/{key}/{key}.log', f'{path}/{key}/{key}_1.log')
                # os.rename(f'{path}/{key}/{key}.chk', f'{path}/{key}/{key}_1.chk')

                
                xyz_list = list()
                match = False
                count = 0
                for line in lines:
                    if 'Standard orientation' in line:
                        match = True
                        count += 1
                        continue
                    elif match and (count >= 5):
                        if len(line.split()) != 6:
                            break
                        xyz = np.array(line.split()[-3:])
                        xyz = xyz.reshape(-1,3)
                        xyz_list.append(xyz) 

                    elif match:
                        count += 1

                coords = np.vstack(xyz_list)
                
                dictionary[key] = [dictionary[key], coords]
                with open(error_list_file, 'w') as file:
                    for key, value in dictionary.items():
                        file.write(f"{key}: {value}\n")



            return func(dictionary, keys_with_value, value_to_find, *args, **kwargs)
        return wrapper
    return decorator


@process_error_("l103")
def err_l103(error_dict, key_list, value_to_find):
    for key in key_list:
        with open(f'{path}/{key}/{key}_1.gjf', 'r') as f1, open(f'{path}/{key}/{key}.gjf', 'w') as f2:
            match=False
            element=[]
            for line in f1:
                if match == False:
                    if 'Title Card Required' in line:
                        f2.write(line)
                        f2.write('\n')
                        match=True
                    elif '#P' in line:
                        f2.write(line.replace('opt', 'opt=cartesian'))
                    else:
                        f2.write(line)
                elif match == True and line.split():
                    if line.split()[0].isalpha():
                        element.append(line.split()[0])
                    else:
                        f2.write(line)

            if error_dict[key][0] == value_to_find:
                for item in zip(element, error_dict[key][1]):
                    x, y, z = map(float, item[1])
                    f2.write(f"{item[0]}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n")
            f2.write('\n')
            
        sub.call(f'rung16 -q amd {path}/{key}/{key}.gjf',shell = True)
'''

def error_reopt(dictionary_path):
    # error_dict = read_dict_from_file(dictionary_path)
    # error_list=
    # keys_with_value = []
    # for key, val in dictionary.items():
    #     if val == value_to_find:
    #         keys_with_value.append(key)

    df=pd.read_csv(error_list_file)
    error_dict={}
    for key in df['Index']:

        counter = 1
        while os.path.exists(f'{gau_path}/{key}/{key}_{counter}.gjf'):
            counter += 1
        
        os.rename(f'{gau_path}/{key}/{key}.gjf', f'{gau_path}/{key}/{key}_{counter}.gjf')
        # try:
        os.rename(f'{gau_path}/{key}/{key}.log', f'{gau_path}/{key}/{key}_{counter}.log')
        os.rename(f'{gau_path}/{key}/{key}.chk', f'{gau_path}/{key}/{key}_{counter}.chk')
        # except:
        #     pass

        out = open(f'{gau_path}/{key}/{key}_{counter}.log', 'r')
        
        # out = open(f'{gau_path}/{key}/{key}_5.log', 'r')
        lines1 = out.readlines()
        
        coords = get_coords(lines1)
        if len(coords) == 0:
            get_gjf_coords(f'{gau_path}/{key}/{key}_{counter}.gjf')
        lines1.reverse()
        lines = []
        for line in lines1:
            if 'Initial Parameters' in line:
                break
            else:
                lines.append(line)
        lines.reverse()
        labels = get_atom_labels(lines)

        error_dict[key] = [str(df.loc[df['Index'] == key, 'error_type']), labels, coords]
        
        counter=1
        with open(f'{gau_path}/{key}/{key}_{counter}.gjf', 'r') as f1, open(f'{gau_path}/{key}/{key}.gjf', 'w') as f2:
            match=False

            for line in f1:
                if 'Title Card Required' in line:
                    f2.write(line)
                    l_i = f1.readline()
                    f2.write(l_i)
                    l_i = f1.readline()
                    f2.write(l_i)
                    break
                elif '#P' in line:
                    pattern = r'opt\S*'
                    matches = re.findall(pattern, line)[0]
                    if error_dict[key][0] == 'l103':
                        f2.write(line.replace(matches, 'opt=cartesian'))
                            
                    elif error_dict[key][0] == 'l9999' or error_dict[key][0] == 'l401':
                        f2.write(line.replace(matches, 'opt=(MaxStep=5, NoTrustUpdate, GDIIS, calcfc)'))
                        
                    elif error_dict[key][0] == 'l701':
                        f2.write(line.replace(matches, 'opt SCF=(noincfock, vshift=300)'))

                else:
                    f2.write(line)

            
            # if error_dict[key][0] == value_to_find:
            
            for items in zip(error_dict[key][1], error_dict[key][2]):
                x, y, z = map(float, items[1])
                f2.write(f"{items[0]}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n")
            f2.write('\n')  

def get_atom_labels(lines):
    atomic_number_to_element = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B',
        6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 16: 'S'
        }

    atoms=[]
    match = False
    count = 0
    for line in lines:
        if 'Standard orientation' in line:
            match = True
            count += 1
            continue
        elif match and (count >= 5):
            if len(line.split()) != 6:
                break
            atoms.append(int(line.split()[1]))
        elif match:
            count += 1
    labels = [atomic_number_to_element.get(number) for number in atoms]
    return labels

def get_coords(lines):
    xyz_list = list()
    match = False
    count = 0
    ene, coords = {}, {}
    num_coord =  0
    for line in lines:
        if 'Standard orientation' in line:
            num_coord+=1
            match = True
            count += 1
            continue
        elif 'SCF Done:' in line:
            ene[num_coord]=float(line.split()[4])
        elif match and (count >= 5):
            if len(line.split()) != 6:
                count = 0
                match = False
                coords[num_coord]=np.vstack(xyz_list)
                xyz_list = list()
                continue
            xyz = np.array(line.split()[-3:], dtype=np.float64)
            xyz = xyz.reshape(-1,3)
            xyz_list.append(xyz)
        elif match:
            count += 1
    if len(ene) != 0:
        min_key = min(ene, key=ene.get)
        coord = coords[min_key]
    else:
        coord=[]
        
    return coord

def get_gjf_coords(file_name):
    with open(file_name, 'r') as f1:
        finding = False
        xyz_list=[]
        for line in f1:
            if finding == True:
                try:
                    xyz = np.array(line.split()[-3:], dtype=np.float64)
                    xyz = xyz.reshape(-1,3)
                    xyz_list.append(xyz)
                except:
                    break
                
            if line[0].isdigit():
                finding = True
               
        coord=np.vstack(xyz_list)   
        
    return coord

def re_rung16_err(error_dict):
    nu = -1
    for i in error_dict:
        nu+=1
        # shutil.copyfile(f'/data/home/zhuyf/dataset_work/database/qm9_database/gau_b3lyp/{i}.gjf',f'/data/home/zhuyf/dataset_work/database/qm9_database/Gaussian_aromatic/{i}/{i}.gjf')
        if nu%5 == 0:
            time.sleep(5)
        os.chdir(f'{gau_path}/{i}')
        sub.call(f'rung16 ./{i}.gjf',shell = True)
        # sub.call(f'mv ./{i}_1.gjf ./{i}.gjf',shell = True)
        # sub.call(f'mv ./{i}_1.chk ./{i}.chk',shell = True)
        # sub.call(f'mv ./{i}_1.log ./{i}.log',shell = True)
        os.chdir(f'{gau_path}')
        
'''
def re_rung16_err_batch(error_dict):
    nu_single, reminder = divmod(len(error_dict), 60)
    for i in range(60):
        os.chdir(f'{gau_path}')
        if i+1 == 60:
            end_nu=(i+1)*nu_single+reminder
        else:
            end_nu=(i+1)*nu_single
        sub.call(f'rung16 -n 2 {{{i*nu_single}..{end_nu}}}/*.gjf',shell = True)
'''



         
if __name__ == '__main__':
    global gau_path
    
    #label = 'aromatic'
    #label = 'bond_but_nonaromatic'
    label = 'double_and_triple'
    work_type = 'excited'

    #gau_path = f'/data/home/zhuyf/dataset_work/database/qm9_database/Gaussian'
    gau_path = f'/data/home/zhuyf/dataset_work/database/qm9_database/Gaussian_gau_excited_wb97xd'
    # input_file="/data/home/zhuyf/dataset_work/database/qm9_database/aromatic.csv"
    input_file = f'/data/home/zhuyf/dataset_work/database/qm9_database/{label}.csv'
    
    error_list_file = f'{gau_path}/error_list_{label}.csv'
    
    input_df=pd.read_csv(input_file)
    if work_type == 'excited':
        gs_opt_file = f'/data/home/zhuyf/dataset_work/database/qm9_database/Gaussian/error_list_{label}.csv'
        error_list= pd.read_csv(gs_opt_error_file)['Index']
        input_df = input_df[~input_df['Index'].isin(error_list)]

    global method, choose_name, choose
    method='b3lyp'
    choose_name = 'gau'

    s = input('1 -- check_log\n2 -- detect_error_type\n3 -- regenerate_gjf\n4 -- re_submit_gau\n')

    if s == '1':
        check_main(input_df, error_list_file)
    elif s == '2':
        error_type(input_df, error_list_file)
    elif s == '3':
        error_reopt(error_list_file)
    elif s == '4':
        error_dict=pd.read_csv(error_list_file)['Index'].tolist()

        if len(error_dict) <= 200:
            re_rung16_err(error_dict)
        else:
            print('\nToo many jobs. Please turn to auto-program!\n')
            # re_rung16_err_batch(error_dict)
            
