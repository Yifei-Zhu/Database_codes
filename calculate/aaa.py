import h5py
import os
import tempfile
from rdkit import Chem
import subprocess as sub

    
with h5py.File('qm9_props.hdf5', 'r') as f:
   # print(f['211'].attrs['SMILES_PYBEL'])
    #print(f['20659']['ground_state']['Etot'][()])
   for i in f:
       print(i)
