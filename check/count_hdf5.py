import h5py
# HDF5 文件路径
hdf5_file_path = 'final_all.hdf5'
i=0
with h5py.File(hdf5_file_path, 'r') as file:
    # 遍历根级别的每个项
    for name in file:
        i+=1     
        # 如果这个项是数据集，增加计数器
        #print(file[name]['ground_state']['Etot'][()])
print(i)
