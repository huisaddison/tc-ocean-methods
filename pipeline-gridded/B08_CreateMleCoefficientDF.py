import argparse
import h5py
import numpy as np
import pandas as pd
import pickle as pkl
from tools import create_ArgoProfileID

parser = argparse.ArgumentParser(description='Create GP coefficient dataframe')
parser.add_argument('--integrated', dest='mode', action='store_const',
                    const='integrated', default='gridded',
                    help='process assuming integrated heat content (default: gridded)')
args = parser.parse_args()

data_dir = 'Data/'
results_dir = 'Results/'
window_size_gp = 5

if args.mode == 'gridded':
    base_fn = (f'{results_dir}localMLESpaceTimeCoefs_Depth_DEPTH_8_20_{window_size_gp}'
               +'_09_2007_2018_AllBasins.mat')

    grid_start = 10
    grid_end = 200
    grid_stride = 10
    n_grid = int(200 / 10)

    depth = 10
    fn_in = base_fn.replace('DEPTH', f'{depth:03d}')
    f = h5py.File(fn_in, 'r')
    FILES = [h5py.File(base_fn.replace('DEPTH', f'{depth:03d}'), 'r')
             for depth in range(grid_start, grid_end+1, grid_stride)]
else:
    n_grid = 1
    FILES = [h5py.File(
        f'{results_dir}localMLESpaceTimeCoefs_8_20_{window_size_gp}'
         '_09_2007_2018_AllBasins.mat')]

# Create Array of ArgoProfileIDs
profile_ids = [
        create_ArgoProfileID(arr_fid[0], arr_cn[0])
        for arr_fid, arr_cn in zip(
            np.array(FILES[0]['FloatIDReg']),
            np.array(FILES[0]['CycleNumberReg']),
        )]
n_prof = len(profile_ids)

# Create mappings from ProfileID to Coefficients
VAR = [
    ('thetas',  'Phi'),
    ('thetat',  'Thetat'),
    ('sigma',   'Sigma'),
]

for var, export_name in VAR:
    mat = np.zeros((n_prof, n_grid))
    for idx in range(n_grid):
        mat[:, idx] = np.array(FILES[idx][var]).flatten()
        dict_ = {
            k: v for k, v in zip(
                profile_ids,
                mat,
        )}
    pkl.dump(dict_, open(
        f'{results_dir}/Mle{export_name}Dict_{window_size_gp}.pkl', 'wb'))

phi_dict = pkl.load(open(
    f'{results_dir}/MlePhiDict_{window_size_gp}.pkl', 'rb'))
theta_t_dict = pkl.load(
        open(f'{results_dir}/MleThetatDict_{window_size_gp}.pkl', 'rb'))
sigma_dict = pkl.load(
        open(f'{results_dir}/MleSigmaDict_{window_size_gp}.pkl', 'rb'))

# Create base DataFrame of unique before_pids
HU = pkl.load(open(f'{data_dir}/HurricaneAdjRawTempDF.pkl', 'rb'))
df = pd.DataFrame({
    'before_pid': HU['before_pid'].unique()
})

# Attach coefficient arrays for before_pids
df['phi'] = df['before_pid'].map(phi_dict)
df['theta_t'] = df['before_pid'].map(theta_t_dict)
df['sigma'] = df['before_pid'].map(sigma_dict)


# Create column to denote whether there is a nan at any row
phi_nan = df['phi'].apply(lambda x: np.sum(np.isnan(x)))
theta_t_nan = df['theta_t'].apply(lambda x: np.sum(np.isnan(x)))
sigma_nan = df['sigma'].apply(lambda x: np.sum(np.isnan(x)))
np.sum(phi_nan>0)
np.sum(theta_t_nan>0)
np.sum(sigma_nan>0)

# Save dataframe
df.set_index('before_pid', inplace=True)
pkl.dump(df, open(f'{results_dir}/MleCoefDF_{window_size_gp}.pkl', 'wb'))

