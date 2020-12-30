import h5py
import numpy as np
import pandas as pd
import pickle as pkl

from tools import (
        conform_lons,
        create_ArgoProfileID,
        matlab_to_datetime,
        matlab_to_datetime_v,
        )

data_dir = '/run/media/addison/BackupsSSD/ocean/Data/'
year_pairs = (
        (2007, 2010),
        (2011, 2014),
        (2015, 2016),
        (2017, 2018),
    )


df_lst = []
for (sy, ey) in year_pairs:
    path = data_dir + f'Argo_data_aggr_{sy}_{ey}.mat'
    f = h5py.File(path)
    lons =              conform_lons(np.array(f['profLongAggr'])).flatten()
    lats =              np.array(f['profLatAggr']).flatten()
    dates =             np.array(f['profJulDayAggr']).flatten()
    datetimes =         np.array([matlab_to_datetime(x) for x in dates])
    float_ids =         np.array(f['profFloatIDAggr']).flatten().astype(int)
    cycle_nums =        np.array(f['profCycleNumberAggr']).flatten()
    profile_ids =       np.array(list(
        create_ArgoProfileID(f, c)
        for f, c in zip(float_ids, cycle_nums)))
    df_lst.append(pd.DataFrame({
        'ProfileID':        profile_ids,
        'Longitude':        lons,
        'Latitude':         lats,
        'ArgoDate':         dates,
        'Timestamp':        datetimes,
        'FloatID':          float_ids,
        'CycleNum':         cycle_nums,
        }))

df = pd.concat(df_lst)
df = df.reset_index(drop=True)
pkl.dump(df, open('Data/ArgoProfileDF.pkl', 'wb'))
