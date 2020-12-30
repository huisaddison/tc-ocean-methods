import matplotlib.pyplot as plt
import numpy as np
import pickle as pkl
import sys
sys.path.append('../implementations/')
from implementation_tools import (
    temp_diff,
)

# Set global variables
output_dir = 'Estimates/'
DEPTH_IDX = 20

data_dir = 'Data/'
results_dir = 'Results/'
window_size_gp=5
stage = 'adj'

sub = 'Hur'

if sub == 'All':
    df = pkl.load(open(f'{data_dir}/HurricaneTempCovDF_Subset.pkl', 'rb'))
else:
    df = pkl.load(open(f'{data_dir}/HurricaneTempCovDF_Subset{sub}.pkl', 'rb'))

y = temp_diff(df, stage, 0)
ymat = np.zeros((1, 20, len(y)))

for depth_idx in range(20):
    ymat[0, depth_idx, :] = temp_diff(df, stage, depth_idx)



fname_small = 'Estimates/B32_LOOCV_Data.pkl'
fname_large = 'Estimates/B35_LOOCV_Data.pkl'

LAMB_SMALL, LOOCV_SMALL = pkl.load(open(fname_small, 'rb'))
LAMB_LARGE, LOOCV_LARGE = pkl.load(open(fname_large, 'rb'))
LAMB = np.concatenate((LAMB_SMALL, LAMB_LARGE[1:]))
LOOCV = np.concatenate((LOOCV_SMALL, LOOCV_LARGE[1:, :, :, :]), axis=0)
n = LOOCV.shape[3]

# 0 - ignored
# 1 - predictions
# 2 - estimated standard deviations
# 3 - diag_H

varmat = np.vstack(df['var']).T
varmat = varmat.reshape(1, *varmat.shape)

loocv_estimates = ((1/varmat) * (
    (LOOCV[:,:,1,:] - ymat) / (1-LOOCV[:,:,3,:]))**2)


error_estimates = loocv_estimates.sum(axis=2)

error_estimates_nonan = error_estimates.copy()
error_estimates_nonan[np.isnan(error_estimates)] = np.inf

lhat = error_estimates_nonan.argmin(axis=0)
LAMB[lhat]

plt.scatter(LAMB[lhat], range(10, 210, 10), marker='x')
plt.plot(LAMB[lhat], range(10, 210, 10))
plt.gca().invert_yaxis()
plt.xscale('log')
plt.xlim(min(LAMB), max(LAMB))
plt.xlabel('$\lambda$')
plt.ylabel('Pressure, dbars ($z$)')
plt.show()


for depth_idx in range(20):
    if depth_idx < 10:
        _ = plt.plot(LAMB, error_estimates[:, depth_idx], label=f'{(depth_idx+1)*10}')
    else:
        _ = plt.plot(LAMB, error_estimates[:, depth_idx], label=f'{(depth_idx+1)*10}', linestyle='dotted')
    _ = plt.scatter(LAMB[lhat[depth_idx]], error_estimates[lhat[depth_idx], depth_idx],
            marker='x')

plt.xscale('log')
plt.xlabel('$\lambda$')
plt.ylabel('LOOCV Error')
plt.axvline(0.65, color='black')
plt.legend(title='Pressure')


for depth_idx in [2, 12]:
    _ = plt.plot(LAMB, error_estimates[:, depth_idx], label=f'{(depth_idx+1)*10}')
    _ = plt.scatter(LAMB[lhat[depth_idx]], error_estimates[lhat[depth_idx], depth_idx],
            marker='x')

plt.xscale('log')
plt.legend()

depth_idx=2
_ = plt.plot(LAMB, error_estimates[:, depth_idx], label=f'{(depth_idx+1)*10}')
_ = plt.scatter(LAMB[lhat[depth_idx]], error_estimates[lhat[depth_idx], depth_idx],
        marker='x')

plt.xscale('log')
plt.legend()




error_mean = LOOCV.mean(axis=2)
error_std = LOOCV.std(axis=2) / np.sqrt(n)
plt.errorbar(LAMB, error_mean[:, idx], yerr=error_std[:, idx])
plt.errorbar(LAMB[:10], error_mean[:10, idx], yerr=error_std[:10, idx])
