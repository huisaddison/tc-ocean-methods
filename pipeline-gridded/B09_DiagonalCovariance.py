'''
author: addison@stat.cmu.edu

Constructs a diagonal covariance matrix in the form of a DataFrame column.
'''
import numpy as np
import pandas as pd
import pickle as pkl

from itertools import product

from tools import covariance_matrix

INTEGRATION_LENGTH = 190
data_dir = 'Data/'
results_dir = 'Results/'
window_size_gp = 5
HU = pkl.load(open(f'{data_dir}/HurricaneAdjRawTempDF.pkl', 'rb'))
n_grid = len(HU.raw_before_temp.iloc[0])
MLE = pkl.load(open(f'{results_dir}/MleCoefDF_{window_size_gp}.pkl', 'rb'))

HU = HU.sort_values(['before_pid', 'after_pid']).reset_index(drop=True)
df_sub = HU.iloc[[0]]
covariance_matrix(df_sub, MLE, 0)

var_mat = np.zeros((HU.shape[0], n_grid))
for idx, depth_idx in product(range(HU.shape[0]), range(n_grid)):
    try:
        var_val = covariance_matrix(HU.iloc[[idx]], MLE, depth_idx)[0][0]
        if n_grid == 1:
            var_val /= (INTEGRATION_LENGTH ** 2)
        var_mat[idx, depth_idx] = var_val
    except TypeError:
        var_mat[idx, depth_idx] = np.nan

HU['var'] = pd.Series((row for row in var_mat))
pkl.dump(HU, open(f'{data_dir}/HurricaneAdjRawTempCovDF.pkl', 'wb'))
