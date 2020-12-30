'''
author: addison@stat.cmu.edu

Construct block diagonal covariance matrices.
'''
import numpy as np
import pandas as pd
import pickle as pkl
import scipy.sparse as sps
import sys
import multiprocessing

from functools import partial
from itertools import product

from tools import covariance_matrix

CPU_COUNT = multiprocessing.cpu_count()

INTEGRATION_LENGTH = 190
data_dir = 'Data/'
df = pkl.load(open(f'{data_dir}/HurricaneAdjRawTempCovDF.pkl', 'rb'))
df.shape
thresh = -4.5
df = df[np.log(df['var'].apply(lambda x: np.min(x))) >= thresh].sort_values([
    'before_pid',
    'after_pid',
])

sub = sys.argv[1]
if sub == 'Hur':
    df = df[df['wind'] >= 64]
elif sub == 'TSTD':
    df = df[df['wind'] <  64]
else:
    print('Assuming all...')
    sub = ''

df.shape
pkl.dump(df, open(f'{data_dir}/HurricaneTempCovDF_Subset{sub}.pkl', 'wb'))


results_dir = 'Results/'
window_size_gp = 5
n_grid = len(df.raw_before_temp.iloc[0])

depth_idx = 1
before_pids = np.sort(df['before_pid'].unique())

MLE = pkl.load(open(f'{results_dir}/MleCoefDF_{window_size_gp}.pkl', 'rb'))

df_list = [df[df['before_pid']==bp] for bp in before_pids]
for depth_idx in range(n_grid):
    print(depth_idx)
    cov = partial(covariance_matrix, df_param=MLE, depth_idx=depth_idx)

    with multiprocessing.Pool(processes=CPU_COUNT) as pool:
        covmat_list = pool.map(cov, df_list)

    with multiprocessing.Pool(processes=CPU_COUNT) as pool:
        premat_list = pool.map(np.linalg.inv, covmat_list)

    C = sps.block_diag(covmat_list)
    P = sps.block_diag(premat_list)
    if n_grid == 1:
        C /= (INTEGRATION_LENGTH ** 2)
        P *= (INTEGRATION_LENGTH ** 2)
    pkl.dump(C, open(f'{results_dir}/BlockCovmat_{window_size_gp}_'
                     f'{(depth_idx+1)*10:03d}{sub}.pkl', 'wb'))
    pkl.dump(P, open(f'{results_dir}/BlockPremat_{window_size_gp}_'
                     f'{(depth_idx+1)*10:03d}{sub}.pkl', 'wb'))
