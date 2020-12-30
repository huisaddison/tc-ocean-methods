import h5py
import numpy as np
import pandas as pd
import pickle as pkl
from tools import create_ArgoProfileID

df = pkl.load(open('./Data/ArgoProfileDF_NoHur.pkl', 'rb'))
NoHurs = set(df['ProfileID'].values)

path = './Data/gridTempProfFiltered_20.mat'
f = h5py.File(path)
float_ids =         np.array(f['profFloatIDAggrSel']).flatten().astype(int)
cycle_nums =        np.array(f['profCycleNumberAggrSel']).flatten()
profile_ids =       list(
    create_ArgoProfileID(f, c)
    for f, c in zip(float_ids, cycle_nums))
mask = np.array([pid in NoHurs for pid in profile_ids], dtype=int)

np.savetxt('./Masks/NoHurMask.csv', mask)
np.savetxt('./Masks/HurricaneProfileMask.csv', 1-mask)


