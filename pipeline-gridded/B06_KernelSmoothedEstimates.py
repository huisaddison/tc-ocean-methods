from itertools import product
import numpy as np
import os
import pandas as pd
import pickle as pkl
import sys
sys.path.append('../implementations/')

from Regressors import KernelSmoother
from implementation_tools import (
    grid,
    temp_diff,
)

# Set global variables
output_dir = 'Estimates/'
STAGE = ('adj', 'raw', 'mf')
h = 0.2 #bandwidth

# Read in dataframe
df_all = pkl.load(open('Data/HurricaneAdjRawTempDF.pkl', 'rb'))
DFS = [
    (df_all, 'combined'),
    (df_all[df_all['wind'] >= 64], 'hurricanes'),
    (df_all[df_all['wind'] <  64], 'tstd'),
]
DEPTH_IDX = len(df_all['adj_before_temp'].iloc[0])



# Create grid for output
bounds_x1 = (-8, +8)
bounds_x2 = (-2, +20)
shape = (100, 400)
test_X = grid(bounds_x1, bounds_x2, *shape).copy()
n_test = test_X.shape[0]
# Save grid 
pkl.dump(test_X, open(f'{output_dir}/test_X.pkl', 'wb'))


for stage, (df, subset) in product(STAGE, DFS):
    fn_out = f'{output_dir}/KernelSmoothedMatx_{h}_{stage}_{subset}.pkl'
    if os.path.exists(fn_out):
        continue
    estimates = np.zeros((n_test, DEPTH_IDX))
    for depth_idx in range(DEPTH_IDX):
        print(stage, depth_idx)
        ks = KernelSmoother(h=h)
        train_X = np.array(df[['standard_signed_angle', 'hurricane_dtd']]).copy()
        if stage == 'mf':
            raw = temp_diff(df, 'raw', depth_idx)
            adj = temp_diff(df, 'adj', depth_idx)
            train_y = (raw - adj).copy()
        else:
            train_y = temp_diff(df, stage, depth_idx)
        ks.fit(train_X, train_y)
        estimates[:, depth_idx] = ks.predict(test_X)
    pkl.dump(estimates,
        open(fn_out, 'wb'))

