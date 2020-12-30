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
sys.path.append('../pipeline-gridded/')

from implementation_tools import (
    grid,
    temp_diff,
)

import plot_config

input_dir = '../pipeline-integrated/Estimates/'
plot_dir = '../pipeline-integrated/Figures/'

fs = 24
df = 2
plt.rcParams['font.family'] = 'Liberation Serif'
#plt.rcParams['font.serif'] = 'Times New Roman'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['mathtext.it'] = 'serif:italic'
plt.rcParams['mathtext.bf'] = 'serif:bold'
plt.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams.update({'font.size': fs})

ypos = -0.3

prefix = '_LOOCV'
lamb = 'AdaptInflate'
block = pkl.load(open(f'{input_dir}/TPS{prefix}_Preds_Block_{lamb}_Hur.pkl', 'rb'))
bmask = pkl.load(open(f'{input_dir}/TPS{prefix}_Mask_Block_{lamb}_Hur.pkl', 'rb'))
ks = pkl.load(open(f'{input_dir}/KernelSmoothedMatx_0.2_adj_hurricanes.pkl', 'rb'))


idx = 0
mat1 = ks[:, idx]
preds_block = block[:, idx]
mat2 = preds_block
mat3 = preds_block.copy()
mat3[~bmask[:, idx].astype(bool)] = np.nan


clim=plot_config.clim
bounds_x1 = (-8, +8)
bounds_x2 = (-2, +20)
shape = (100, 400)
minima, maxima = -clim, +clim

norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=False)
cmap = cm.bwr
cmap.set_bad(color='gray')
mapper = cm.ScalarMappable(norm=norm, cmap=cmap)

fig = plt.figure(figsize=(20, 12))
gs= gridspec.GridSpec(1, 3, figure=fig,
    width_ratios=[1.0, 1.0, 1.0666])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])

ax1.imshow(mat1.reshape(*(shape[::-1])),
        origin='lower',
        cmap=cmap, norm=norm,
        extent=(*bounds_x1, *bounds_x2),
        )
ax1.invert_yaxis()
ax1.axhline(0, color='k', linewidth=0.5)
ax1.axvline(0, color='k', linewidth=0.5)
ax1.set_xticks([-5, 0, +5])
ax1.set_yticks([0, 5, 10, 15, 20])
ax1.set_ylabel(r'Days since TC passage ($\tau$)',
        fontsize=fs)
ax1.set_title(r'(a) Seasonally adjusted',
        y=ypos,
        fontsize=fs+df)


ax2.imshow(mat2.reshape(*(shape[::-1])),
        origin='lower',
        cmap=cmap, norm=norm,
        extent=(*bounds_x1, *bounds_x2),
        )
ax2.invert_yaxis()
ax2.axhline(0, color='k', linewidth=0.5)
ax2.axvline(0, color='k', linewidth=0.5)
ax2.set_xticks([-5, 0, +5])
ax2.set_yticks([])
ax2.set_xlabel(r'Cross-track angle, days ($d$)',
        fontsize=fs)
ax2.set_title(r'(b) Block covariance',
        y=ypos,
        fontsize=fs+df)

ax3.imshow(mat3.reshape(*(shape[::-1])),
        origin='lower',
        cmap=cmap, norm=norm,
        extent=(*bounds_x1, *bounds_x2),
        )
ax3.invert_yaxis()
ax3.axhline(0, color='k', linewidth=0.5)
ax3.axvline(0, color='k', linewidth=0.5)
ax3.set_xticks([-5, 0, +5])
#ax3.xaxis.set_ticks_position('top') 
ax3.set_title(r'(c) Pointwise $\alpha=0.5$ test',
        y=ypos,
        fontsize=fs+df)

divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(mapper, ax=ax3, cax=cax)
cbar.set_label(r'Temperature difference ($^\circ$C)',
        fontsize=fs)
ax3.set_yticks([])


plt.savefig(f'{plot_dir}/Integrated_ThreePanel.pdf',
        bbox_inches='tight', pad_inches=0)

