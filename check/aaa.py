import h5py

    
with h5py.File('final_all.hdf5', 'r') as f:
    print(f['Ba022862193'].attrs['SMILES_PYBEL'])
    #print(f['20659']['ground_state']['Etot'][()])
   #for i in f:
   #   print(i)
