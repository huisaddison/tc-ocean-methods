from tools import create_ArgoProfileID
import numpy as np
import pandas as pd
import pickle as pkl

def attach_temps(df, dd, prefix=''):
    df[f'{prefix}before_temp'] = df['before_pid'].map(dd)
    df[f'{prefix}after_temp']  = df['after_pid'].map(dd)
    return df

def subsample(df: pd.DataFrame, n_days: int = 3) -> pd.DataFrame:
    '''
    Enforce that subsequent floats in the lineage of a before
    profile are separated by at least n_days.

    Parameters
    ----------
    df: pd.DataFrame
        Dataframe to have rows filtered.  All profiles must have the same
        before float.
    n_days: int
        Number of days subsequent floats must be separated by.

    Returns
    -------
        pd.DataFrame
    '''
    if len(df['before_pid'].unique()) != 1:
        raise ValueError('Must have exactly one before_pid')
    n = df.shape[0]
    keep = np.zeros(n).astype(bool)
    for idx in range(n):
        if (keep.sum() == 0):
            prev_t = df.iloc[0]['before_t']
        else:
            prev_t = df.iloc[max(keep*np.arange(n))]['after_t']
        keep[idx] = (df.iloc[idx]['after_t'] - prev_t
                     >= pd.Timedelta(days=n_days))
    return df[keep]



# Read in dictionaries of values
d_hur = pkl.load(open('./Data/HurricaneResDict.pkl', 'rb'))
d_raw = pkl.load(open('./Data/HurricaneProfDict.pkl', 'rb'))

# Read in database of hurricane profile pairs
HU = pkl.load(open('./Data/AllBasin_PairDF.pkl', 'rb'))
HU = HU.sort_values(['before_pid', 'after_t'])
print(HU.shape)
## 19668 rows

# Filter out the repeats
df_list = []
before_pids = HU['before_pid'].unique()
for bp in before_pids:
    df_list.append(
            subsample(HU[HU['before_pid'] == bp].sort_values('after_t')))

HU = pd.concat(df_list)
print(HU.shape)
## 16896


# Attach temperature information
HU = attach_temps(HU, d_hur, 'adj_')
HU = attach_temps(HU, d_raw, 'raw_')

# Add standard_signed_angle where SH is flipped
signs = HU['sign'].values
is_SH = HU['HurricaneID'].apply(lambda x: x[:2] == 'SH').values
mask_SH = (~is_SH) * 2 - 1
HU['standard_signed_angle'] = HU['signed_angle'] * mask_SH

print(HU.shape)
# 16896

HU = HU[~HU['adj_before_temp'].isna()]
HU = HU[~HU['adj_after_temp'].isna()]
HU = HU[~HU['raw_before_temp'].isna()]
HU = HU[~HU['raw_after_temp'].isna()]
print(HU.shape)
# 16025

HU['hurricane_dtd'] = HU['hurricane_dt'] /  pd.to_timedelta(1, unit='D')
HU['profile_dtd'] = HU['profile_dt'] /  pd.to_timedelta(1, unit='D')

HU['hurricane_id'] = HU['HurricaneID']

# Hack to handle integrated rescaling
if len(HU['adj_before_temp'].iloc[0]) == 1:
    INTEGRATION_LENGTH = 190
    print('Value arrays of length 1; assuming integration length of '
            f'{INTEGRATION_LENGTH} and rescaling...')
    HU['adj_before_temp'] = HU['adj_before_temp'] / INTEGRATION_LENGTH
    HU['adj_after_temp'] = HU['adj_after_temp'] / INTEGRATION_LENGTH
    HU['raw_before_temp'] = HU['raw_before_temp'] / INTEGRATION_LENGTH
    HU['raw_after_temp'] = HU['raw_after_temp'] / INTEGRATION_LENGTH

pkl.dump(HU, open('Data/HurricaneAdjRawTempDF.pkl', 'wb'))
