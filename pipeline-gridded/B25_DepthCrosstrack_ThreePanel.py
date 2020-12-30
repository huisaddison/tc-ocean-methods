import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import pickle as pkl
import os
import sys
sys.path.append('../implementations/')
import plot_config

fs = 24
df = 2
plt.rcParams['font.family'] = 'Liberation Serif'
#plt.rcParams['font.serif'] = 'Times New Roman'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['mathtext.it'] = 'serif:italic'
plt.rcParams['mathtext.bf'] = 'serif:bold'
plt.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams.update({'font.size': fs})

from implementation_tools import (
    grid,
    temp_diff,
)

input_dir = '../pipeline-gridded/Estimates/'
plot_dir = '../pipeline-gridded/Figures_TPS_ThreePanel/'


prefix = '_LOOCV'
lamb = 'AdaptInflate'

tps_est_dir = 'Estimates/'

BASIS = pkl.load(open(f'{tps_est_dir}/TPS{prefix}_Basis_Block_{lamb}_Hur.pkl',
    'rb'))
THETA = pkl.load(open(f'{tps_est_dir}/TPS{prefix}_Theta_Block_{lamb}_Hur.pkl',
    'rb'))
COV = pkl.load(open(f'{tps_est_dir}/TPS{prefix}_CovTheta_Block_{lamb}_Hur.pkl',
    'rb'))

def depthxc_estimates(d1, d2,
        BASIS=BASIS, THETA=THETA, COV=COV):
    n_params = THETA.shape[0]

    bounds_x1 = (-8, +8)
    bounds_x2 = (-2, +20)
    shape = (100, 400)
    test_X = grid(bounds_x1, bounds_x2, *shape)

    xc_filt = ((test_X[:, 1] >= d1)
              *(test_X[:, 1] < d2))

    D = np.sort(np.unique(test_X[:, 0]))
    n_d = len(D)
    DEPTH_IDX = 20

    predmat = np.zeros((DEPTH_IDX, n_d))
    maskmat = np.zeros((DEPTH_IDX, n_d))

    for tidx, d in enumerate(D):
        idxs = xc_filt * (test_X[:, 0] == d)
        x_pts = test_X[idxs, 1]  # 1 x 6
        y_pts = BASIS[idxs, :]   # 6 x 394
        # trap integration for each column
        int_basis = np.zeros(n_params)
        var = np.zeros(DEPTH_IDX)
        for pidx in range(n_params):
            int_basis[pidx] = np.trapz(y_pts[:, pidx], x_pts)
        for depth_idx in range(DEPTH_IDX):
            var[depth_idx] = np.linalg.multi_dot((int_basis,
                                COV[:, :, depth_idx],
                                int_basis))
        inprod = np.dot(int_basis, THETA)
        predmat[:, tidx] = inprod / (max(x_pts) - min(x_pts))
        sd = np.sqrt(var)
        maskmat[:, tidx] = (inprod > 2 * sd) | (inprod < -2 * sd)
    predmat_small = predmat
    maskmat_small = maskmat
    predmat = np.zeros((DEPTH_IDX, n_d*4))
    maskmat = np.zeros((DEPTH_IDX, n_d*4))
    xnew = np.linspace(-8, +8, n_d*4)
    for depth_idx in range(DEPTH_IDX):
        predmat[depth_idx, :] = np.interp(xnew, D, predmat_small[depth_idx, :])
        maskmat[depth_idx, :] = np.repeat(maskmat_small[depth_idx, :], 4)
    mat = predmat.copy()
    mat[~maskmat.astype(bool)] = np.nan
    return mat


try:
    mat1, mat2 = pkl.load(open('cache/tmp_depthxc.pkl', 'rb'))
except FileNotFoundError:
    mat1 = depthxc_estimates(0, 3)
    mat2 = depthxc_estimates(3, 20)
    pkl.dump((mat1, mat2), open('cache/tmp_depthxc.pkl', 'wb'))

bounds_x1 = (-8, +8)
bounds_x2 = (5, 205)

clim=plot_config.clim
shape = (100, 400)
minima, maxima = -clim, +clim

norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=False)
cmap = cm.bwr
cmap.set_bad(color='gray')
mapper = cm.ScalarMappable(norm=norm, cmap=cmap)

fig = plt.figure(figsize=(20, 6))
gs= gridspec.GridSpec(1, 2, figure=fig,
    width_ratios=[1.0, 1.0006])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

ax1.imshow(mat1,
        origin='lower',
        cmap=cmap, norm=norm,
        extent=(*bounds_x1, *bounds_x2),
        aspect='auto',
        )
ax1.invert_yaxis()
ax1.axvline(0, color='k', linewidth=0.5)
ax1.set_xticks([-5, 0, 5])
ax1.set_yticks([10, 50, 100, 150, 200])
ax1.set_ylabel(r'Pressure, dbars ($z$)',
        fontsize=fs)
'''
ax1.set_xlabel(r'Cross-track angle, degrees ($d$)',
        fontsize=fs)
'''
ypos = -0.35
ax1.set_title(r'(a) $\tau \in [0, +3)$ days',
        y=ypos,
        fontsize=fs+df)


ax2.imshow(mat2,
        origin='lower',
        cmap=cmap, norm=norm,
        extent=(*bounds_x1, *bounds_x2),
        aspect='auto',
        )
ax2.invert_yaxis()
ax2.axvline(0, color='k', linewidth=0.5)
ax2.set_xticks([-5, 0, 5])
ax2.set_yticks([])
#ax2.xaxis.set_ticks_position('top') 
'''
ax2.set_xlabel(r'Cross-track angle, degrees ($d$)',
        fontsize=fs)
'''
ax2.set_title(r'(b) $\tau \in [+3, +20)$ days',
        y=ypos,
        fontsize=fs+df)

divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(mapper, ax=ax2, cax=cax)
cbar.set_label(r'Temperature difference ($^\circ$C)',
        fontsize=fs)

fig.text(0.5, -0.03, r'Cross-track angle, degrees ($d$)',
        ha='center',
        fontsize=fs)

plt.savefig(f'{plot_dir}/TPS_ThreePanel_DepthXC.pdf',
        bbox_inches='tight', pad_inches=0)
