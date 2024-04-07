import os
import re
import sys
import numpy as np
import pandas as pd
import pickle
import periodictable
from openbabel import pybel
from typing import Any, Iterator, List, Optional, Tuple, Union, cast, IO
from collections import defaultdict
import ast

# from scipy.sparse import lil_matrix
# from multiprocessing import Pool
# def porcess_ECFP(fp_list):
#     num_processes=8
#     num_strings = len(fp_list)
#     string_length = len(fp_list[0])
#     binary_matrix = lil_matrix((num_strings, string_length), dtype=int)
#     with Pool(processes=num_processes) as pool:
#         results = pool.starmap(process_string, [(fp_list[i], i) for i in range(num_strings)])
#     for index, binary_vector in results:
#         binary_matrix[index, :] = binary_vector
#     return binary_matrix.toarray()

# def process_string(string, index):
#     binary_vector = [int(bit) for bit in string]
#     return index, binary_vector

def readECFP_single(ind):
    loaded_data = np.load(f'{gau_path}/{ind}/ecfp_data.bin')
    return loaded_data

def generateECFP4Matrix(main_file):
    fp_list= []
    df_ind = pd.read_csv(main_file)['Index']
    for i, ind in enumerate(df_ind):
        if ind != 129152:
            ecfp=readECFP_single(ind)
            ecfp = np.array([int(bit) for bit in str(ecfp)], dtype=np.uint8)
            fp_list.append(ecfp)
    fp_matrix = np.vstack(fp_list)
    fp_matrix = keepNonConstantCol(fp_matrix)
    #  fp_matrix = porcess_ECFP(fp_list)
    np.save(f'{descriptor_path}/ECFP4_{np.shape(fp_matrix)[0]}samples_{np.shape(fp_matrix)[1]}features_.npy',fp_matrix)

def pad_array(x: np.ndarray,shape: Union[Tuple, int],
    fill: float = 0.0, both: bool = False) -> np.ndarray:
    """
    Pad an array with a fill value.

    Parameters
    ----------
    x: np.ndarray
        A numpy array.
    shape: Tuple or int
        Desired shape. If int, all dimensions are padded to that size.
    fill: float, optional (default 0.0)
        The padded value.
    both: bool, optional (default False)
        If True, split the padding on both sides of each axis. If False,
        padding is applied to the end of each axis.

    Returns
    -------
    np.ndarray
        A padded numpy array
    """
    x = np.asarray(x)
    if not isinstance(shape, tuple):
        shape = tuple(shape for _ in range(x.ndim))
    pad = []
    for i in range(x.ndim):
        diff = shape[i] - x.shape[i]
        assert diff >= 0
        if both:
            a, b = divmod(diff, 2)
            b += a
            pad.append((a, b))
        else:
            pad.append((0, diff))
    pad = tuple(pad)  # type: ignore
    x = np.pad(x, pad, mode='constant', constant_values=fill)
    return x


def randomize_coulomb_matrix(m: np.ndarray, n_samples) -> List[np.ndarray]:
    """Randomize a Coulomb matrix as decribed in [1]_:

    1. Compute row norms for M in a vector row_norms.
    2. Sample a zero-mean unit-variance noise vector e with dimension
        equal to row_norms.
    3. Permute the rows and columns of M with the permutation that
        sorts row_norms + e.

    Parameters
    ----------
    m: np.ndarray
        Coulomb matrix.

    Returns
    -------
    List[np.ndarray]
        List of the random coulomb matrix

    References
    ----------
    .. [1] Montavon et al., New Journal of Physics, 15, (2013), 095003

    """
    rval = []
    row_norms = np.asarray([np.linalg.norm(row) for row in m], dtype=float)
    rng = np.random.RandomState(14)
    for i in range(n_samples):
        e = rng.normal(size=row_norms.size)
        p = np.argsort(row_norms + e)
        new = m[p][:, p]  # permute rows first, then columns
        rval.append(new)
    return rval

def generateOriginalCMMatrices(main_file):
    max_nu = 10+(2*10+2)
    random_option = True

    fp_list = []
    df_ind = pd.read_csv(main_file)['Index']
    for i, ind in enumerate(df_ind):
        if ind != 129152:
            path_i = f'{gau_path}/{ind}'
            cm_file = f'{path_i}/cm.dat'
            cm = np.loadtxt(cm_file)
            if random_option == True:
                nu_cm =5
                for random_m in randomize_coulomb_matrix(cm,nu_cm):
                    random_m = pad_array(random_m, max_nu)
                    fp_list.append(random_m)
            else:
                cm_m = pad_array(cm, max_nu)
                fp_list.append(cm_m)
    return np.asarray(fp_list)

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def apply_sigmoid(row):
    #hyperparameters
    Theta=1
    N_dim=3
    ######
    return np.column_stack([sigmoid(row + i * Theta) for i in range(-int(np.floor(N_dim/2)), int(np.ceil(N_dim/2)))])



def binarizationVecs(cm_vecs):
    return np.apply_along_axis(apply_sigmoid, axis=1, arr=cm_vecs)

def reshapeTensor(m):
    return m.reshape(np.shape(m)[0], np.shape(m)[1]*np.shape(m)[2])

def generateCMMatrix_main(main_file):
    cm_matrices = generateOriginalCMMatrices(main_file)
    print(np.shape(cm_matrices))
    cm_vecs = reshapeTensor(cm_matrices)
    print(np.shape(cm_vecs))
    bi_cm_vecs = binarizationVecs(cm_vecs)
    print(np.shape(bi_cm_vecs))
    bi_cm_vecs = reshapeTensor(bi_cm_vecs)
    print(np.shape(bi_cm_vecs))
    bi_cm_vecs = keepNonConstantCol(bi_cm_vecs)
    print(np.shape(bi_cm_vecs))
    np.save(f'{descriptor_path}/CM_{np.shape(bi_cm_vecs)[0]}samples_{np.shape(bi_cm_vecs)[1]}features.npy',bi_cm_vecs)

def keepNonConstantCol(m):
    non_constant_columns = np.where(~np.all(m == m[0, :], axis=0))[0]
    return m[:, non_constant_columns]


def extractOnePropsFromAll(ind,props):
    file_GS = f'{gau_path}/{ind}/GS_props.dat'
    with open(file_GS,'r') as f1:
        for line in f1:
            if 'PROP:' in line:
                prop_i=''
                key_prop = line.split(':')[-1].strip()
                continue
            if key_prop == "labels":
                continue
            # if ']]' in line:
            #     finding=False
            new_line = line.replace('[','').replace(']','').strip()
            prop_i+= ' '+new_line

            if ']]' in line:
                props[key_prop].append(prop_i.split())
                finding = False

def save_dict(dict_i):
    with open(f'{descriptor_path}/all_props_dict.pkl', 'wb') as pickle_file:
        pickle.dump(dict_i, pickle_file)
def read_dict(pkl_file):
    with open(pkl_file, 'rb') as pickle_file:
        props = pickle.load(pickle_file)

def generateGSTaskVec(main_file):
    # all_GS_properties = ['labels','coords','Etot','e_homo_lumo','polarizability','dipole','quadrupole','zpve','rot_constants','elec_spatial_ext','forces','e_thermal','freqs','mulliken','cv']
    props= defaultdict(list)
    df_ind = pd.read_csv(main_file)['Index']
    for ind in df_ind:
        if ind != 129152:
            extractOnePropsFromAll(ind,props)
    save_dict(props)

    for k in props.keys():
        float_data = [[float(value) for value in sublist] for sublist in props[k]]
        max_length = max(len(sublist) for sublist in float_data)
        filled_data = [sublist + [0] * (max_length - len(sublist)) for sublist in float_data]
        matrix = np.array(filled_data)
        np.save(f'{descriptor_path}/{k}_{np.shape(matrix)[0]}samples_{np.shape(matrix)[1]}features.npy',matrix)




def main():
    global code_path, main_path, gau_path, descriptor_path

    code_path = '/data/home/zhuyf/dataset_work/database/analysis'
    main_path = '/data/home/zhuyf/dataset_work/database/qm9_database'
    gau_path = f'{main_path}/Gaussian'

    descriptor_path = f'{code_path}/descriptor'
    os.makedirs(descriptor_path, exist_ok=True)

    main_file = f'{main_path}/double_and_triple.csv'
    # main_file = f'{main_path}/a.csv'

    generateECFP4Matrix(main_file)
    generateCMMatrix_main(main_file)
    generateGSTaskVec(main_file)

if __name__ == '__main__':
    main()
    # main_test_smi_inchi()
    print('Done!\n')


