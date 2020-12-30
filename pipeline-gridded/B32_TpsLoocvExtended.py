from itertools import product
import argparse
import numpy as np
import os
import pandas as pd
import pickle as pkl
import sys
sys.path.append('../implementations/')

from Regressors import ThinPlateSpline
from implementation_tools import (
    grid,
    temp_diff,
)

parser = argparse.ArgumentParser(
        description='Plot thin plate spline estimates (three panel)')
parser.add_argument('--integrated', dest='mode', action='store_const',
                    const='integrated', default='gridded',
                    help='process assuming integrated heat content (default: gridded)')
args = parser.parse_args()

# Set global variables
output_dir = 'Estimates/'
if args.mode == 'gridded':
    DEPTH_IDX = 20
else:
    DEPTH_IDX = 1

data_dir = 'Data/'
results_dir = 'Results/'
window_size_gp=5
stage = 'adj'

sub = 'Hur'

if sub == 'All':
    df = pkl.load(open(f'{data_dir}/HurricaneTempCovDF_Subset.pkl', 'rb'))
else:
    df = pkl.load(open(f'{data_dir}/HurricaneTempCovDF_Subset{sub}.pkl', 'rb'))

X = np.array(df[['standard_signed_angle', 'hurricane_dtd']]).copy()

bounds_x1 = (-8, +8)
bounds_x2 = (-2, +20)
train_knots = grid(bounds_x1, bounds_x2, 33, 45)
n_param = train_knots.shape[0] + 3
shape = (100, 400)
test_X = grid(bounds_x1, bounds_x2, *shape)
n_test = test_X.shape[0]



LAMB_ = np.linspace(0.05, 10.0, 200) 
LAMB = np.zeros(len(LAMB_)+1)
LAMB[0] = 0.01
LAMB[1:] = LAMB_

y = temp_diff(df, stage, 0)
LOOCV = np.zeros((len(LAMB), DEPTH_IDX, 4, len(y)))

for depth_idx in range(DEPTH_IDX):
    print(depth_idx)
    y = temp_diff(df, stage, depth_idx)
    S = df['var'].apply(lambda x: x[depth_idx]).values
    W = 1 / S
    if sub == 'All':
        block_S = pkl.load(open(f'{results_dir}/BlockCovmat_{window_size_gp}_'
                                f'{(depth_idx+1)*10:03d}.pkl', 'rb'))
        block_W = pkl.load(open(f'{results_dir}/BlockPremat_{window_size_gp}_'
                                f'{(depth_idx+1)*10:03d}.pkl', 'rb'))
    else:
        block_S = pkl.load(open(f'{results_dir}/BlockCovmat_{window_size_gp}_'
                                f'{(depth_idx+1)*10:03d}{sub}.pkl', 'rb'))
        block_W = pkl.load(open(f'{results_dir}/BlockPremat_{window_size_gp}_'
                                f'{(depth_idx+1)*10:03d}{sub}.pkl', 'rb'))
     
    for idx, lamb in enumerate(LAMB):
        tps_block = ThinPlateSpline(lamb=lamb, knots=train_knots)
        tps_block.fit(X, y, W=block_W)
        ret = tps_block.predict(X, sd='diag', S=block_S, k=2, diag_H=True)
         
        var_inv = (1/ret[1])**2
        norm_var_inv = var_inv / var_inv.sum()
        resid = (y - ret[0])
        diag_H = ret[3]
         
        LOOCV[idx, depth_idx, 0, :] = norm_var_inv * (
                (resid / (1 - diag_H)) ** 2)
        LOOCV[idx, depth_idx, 1, :] = ret[0]
        LOOCV[idx, depth_idx, 2, :] = ret[1]
        LOOCV[idx, depth_idx, 3, :] = ret[3]

loocv_data = (LAMB, LOOCV)
pkl.dump(loocv_data, open('Estimates/B32_LOOCV_Data.pkl', 'wb'))

