import h5py

hdf5_file = 'final_all.hdf5'

with h5py.File(hdf5_file, 'r') as file:
    for dataset_name in file:
        dataset = file[dataset_name]
        smiles = dataset.attrs['SMILES']
        if smiles not in uni_list:
            uni_list.append(smiles)
        else:
            # dup_list.append(smiles)
            dup_name_list.append(dataset_name)

return dup_name_list


datasets = {
        'pubchemqc9' : {
            'prefix': 'Ba',

            'main_path': '/share/home/limg/analysis_lmg',
            'gs_dir': '/share/home/limg/pubchemqc/Gaussian_1_10_CNOF',
            'es_dir': '/data/home/zhuyf/pubchem_1_10_CNOF/Gaussian_1_10_CNOF_gau_excited_wb97xd',

            'csv_file': '9_CNOF_g_s_2.csv',
            'gs_error': 'pub_10_CNOF_error_list.csv',
            'es_error':'error_list_.csv',

            'hdf5':'9_CNOF_pub_molecules.hdf5',
        },
        'pubchemqc10' : {
            'prefix': 'Bb',
            'main_path': '/share/home/limg/analysis_lmg',
            'gs_dir': '/share/home/limg/pubchemqc/Gaussian_1_10_CNOF',
            'es_dir': '/data/home/zhuyf/pubchem_1_10_CNOF/Gaussian_1_10_CNOF_gau_excited_wb97xd',

            'csv_file': '10_CNOF_g_s_2.csv',
            'gs_error':'pub_10_CNOF_error_list.csv',
            'es_error':'error_list_.csv',

            'hdf5':'10_CNOF_pub_molecules.hdf5',
        },
        'qm9' : {
            'prefix': 'Aa',

            'main_path': '/data/home/zhuyf/dataset_work/database/qm9_database',
            'gs_dir': '/data/home/zhuyf/dataset_work/database/qm9_database/Gaussian',
            'es_dir': '/data/home/zhuyf/dataset_work/database/qm9_database/Gaussian_gau_excited_wb97xd',

            'csv_file': 'double_and_triple.csv',
            'gs_error':'error_list_double_and_triple.csv',
            'es_error': 'error_list_double_and_triple.csv',

            'hdf5':'qm9_props.hdf5',
        },
        'gdb11' : {
            'prefix': 'Ab',

            'main_path': '/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation',
            'gs_dir': '/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation/calculation_/gdb11_Gaussian_',
            'es_dir': '/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation/calculation_excited/gdb11_Gaussian_excited',

            'csv_file': 'new_mbkmeans_final_sample_1w.csv',
            'gs_error':'error_list_.csv',
            'es_error': 'error_list_excited.csv',

            'hdf5':'gdb11_props.hdf5',
        },
    }

