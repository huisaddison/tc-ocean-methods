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
lamb_choice = 'Mean'
lamb_choice = 'Median'
lamb_choice = 'Adapt'
lamb_choice = 10
lamb_choice = 'AdaptSignal'
lamb_choice = 'AdaptSignal2'
lamb_choice = 'AdaptSignalInflate'
lamb_choice = 'AdaptInflate'

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

y = temp_diff(df, stage, 0)
ymat = np.zeros((1, DEPTH_IDX, len(y)))

for depth_idx in range(DEPTH_IDX):
    ymat[0, depth_idx, :] = temp_diff(df, stage, depth_idx)

fname_small = 'Estimates/B32_LOOCV_Data.pkl'
fname_large = 'Estimates/B35_LOOCV_Data.pkl'

LAMB_SMALL, LOOCV_SMALL = pkl.load(open(fname_small, 'rb'))
LAMB_LARGE, LOOCV_LARGE = pkl.load(open(fname_large, 'rb'))
LAMB = np.concatenate((LAMB_SMALL, LAMB_LARGE[1:]))
LOOCV = np.concatenate((LOOCV_SMALL, LOOCV_LARGE[1:, :, :, :]), axis=0)
varmat = np.vstack(df['var']).T
varmat = varmat.reshape(1, *varmat.shape)
loocv_estimates = ((1/varmat) * (
    (LOOCV[:,:,1,:] - ymat) / (1-LOOCV[:,:,3,:]))**2)
error_estimates = loocv_estimates.sum(axis=2)
error_estimates_nonan = error_estimates.copy()
error_estimates_nonan[np.isnan(error_estimates)] = np.inf
lhat = error_estimates_nonan.argmin(axis=0)
argmin_lamb = LAMB[lhat]


inflation_factor=1.01
if lamb_choice == 'Adapt':
    lamb_seq = argmin_lamb
elif lamb_choice == 'Mean':
    lamb_seq = [np.mean(argmin_lamb) for _ in range(DEPTH_IDX)]
elif lamb_choice == 'Median':
    lamb_seq = [np.median(argmin_lamb) for _ in range(DEPTH_IDX)]
elif lamb_choice == 'AdaptSignal':
    signal_region = (df.angle <= 3).values * (0 <= df.hurricane_dtd).values * (df.hurricane_dtd <= 12).values
    error_estimates_signal = loocv_estimates[:,:,signal_region].sum(axis=2)
    error_estimates_nonsignal = loocv_estimates[:,:,~signal_region].sum(axis=2)
    lhat_signal = error_estimates_signal.argmin(axis=0)
    lhat_nonsignal = error_estimates_nonsignal.argmin(axis=0)
    lamb_seq = LAMB[lhat_signal]
elif lamb_choice == 'AdaptSignal2':
    signal_region = (df.angle <= 2).values * (0 <= df.hurricane_dtd).values * (df.hurricane_dtd <= 12).values
    error_estimates_signal = loocv_estimates[:,:,signal_region].sum(axis=2)
    error_estimates_nonsignal = loocv_estimates[:,:,~signal_region].sum(axis=2)
    lhat_signal = error_estimates_signal.argmin(axis=0)
    lhat_nonsignal = error_estimates_nonsignal.argmin(axis=0)
    lamb_seq = LAMB[lhat_signal]
elif lamb_choice == 'AdaptSignalInflate':
    signal_region = ((df.angle <= 2).values * (0 <= df.hurricane_dtd).values
            * (df.hurricane_dtd <= 12).values)
    error_estimates_signal = loocv_estimates[:,:,signal_region].sum(axis=2)
    min_error = error_estimates_signal.min(axis=0)
    lhat_inflated_signal = np.argmax(error_estimates_signal < min_error *
	    inflation_factor, axis=0)
    lamb_seq = LAMB[lhat_inflated_signal]
elif lamb_choice == 'AdaptInflate':
    min_error = error_estimates.min(axis=0)
    inflation_factor = 1.01
    # find the smallest lambda that has error smaller
    lhat_inflated_all = np.argmax(error_estimates < min_error * inflation_factor, axis=0)
    lamb_seq = LAMB[lhat_inflated_all]
else:
    lamb_seq = [lamb_choice for _ in range(DEPTH_IDX)]

PREDS_NOREW = np.zeros((n_test, DEPTH_IDX))
PREDS_BLOCK = np.zeros((n_test, DEPTH_IDX))
STDEV_BLOCK = np.zeros((n_test, DEPTH_IDX))
MASK_BLOCK = np.zeros((n_test, DEPTH_IDX))

THETA_BLOCK = np.zeros((n_param, DEPTH_IDX))
BASIS_BLOCK = np.zeros((n_test, n_param)) # Invariant to depth
COV_THETA_BLOCK = np.zeros((n_param, n_param, DEPTH_IDX))

for depth_idx, lamb in zip(range(DEPTH_IDX), lamb_seq):
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

    tps_norew = ThinPlateSpline(lamb=lamb, knots=train_knots)
    tps_block = ThinPlateSpline(lamb=lamb, knots=train_knots)
    tps_norew.fit(X, y)
    tps_block.fit(X, y, W=block_W)
    ret0 = tps_norew.predict(test_X)
    ret2 = tps_block.predict(test_X, sd='diag', S=block_S, k=2)
    PREDS_NOREW[:, depth_idx] = ret0
    PREDS_BLOCK[:, depth_idx] = ret2[0]
    STDEV_BLOCK[:, depth_idx] = ret2[1]
    MASK_BLOCK[ :, depth_idx] = ret2[2]
    THETA_BLOCK[:, depth_idx] = tps_block.theta
    BASIS_BLOCK[:, :] =      tps_block._test_basis
    COV_THETA_BLOCK[:, :, depth_idx] =  tps_block._cov_theta

pkl.dump(PREDS_NOREW,
        open(f'{output_dir}/TPS_LOOCV_Preds_NoRew_{lamb_choice}_{sub}.pkl', 'wb'))
pkl.dump(PREDS_BLOCK,
        open(f'{output_dir}/TPS_LOOCV_Preds_Block_{lamb_choice}_{sub}.pkl', 'wb'))
pkl.dump(STDEV_BLOCK,
        open(f'{output_dir}/TPS_LOOCV_Stdev_Block_{lamb_choice}_{sub}.pkl', 'wb'))
pkl.dump(MASK_BLOCK,
        open(f'{output_dir}/TPS_LOOCV_Mask_Block_{lamb_choice}_{sub}.pkl', 'wb'))

pkl.dump(THETA_BLOCK,
        open(f'{output_dir}/TPS_LOOCV_Theta_Block_{lamb_choice}_{sub}.pkl', 'wb'))
pkl.dump(BASIS_BLOCK,
        open(f'{output_dir}/TPS_LOOCV_Basis_Block_{lamb_choice}_{sub}.pkl', 'wb'))
pkl.dump(COV_THETA_BLOCK,
        open(f'{output_dir}/TPS_LOOCV_CovTheta_Block_{lamb_choice}_{sub}.pkl', 'wb'))

