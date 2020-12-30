import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NULL = object()

def grid(x1, x2, res1, res2):
    x1_grid = np.linspace(*x1, num=res1)
    x2_grid = np.linspace(*x2, num=res2)
    X1, X2 = np.meshgrid(x1_grid, x2_grid)
    newX = np.vstack([X1.flatten(), X2.flatten()]).T
    return newX

def plot_2d(preds, shape, mask=NULL, clim=2, cbar=True):
    bounds_x1 = (-8, +8)
    bounds_x2 = (-2, +20)
    minima, maxima = -clim, +clim
    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=False)
    cmap = cm.bwr
    cmap.set_bad(color='gray')
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    if mask is not NULL:
        preds = preds.copy()
        preds[~mask] = np.nan
    plt.imshow(preds.reshape(*(shape[::-1])),
            origin='lower',
            cmap=cmap, norm=norm,
            extent=(*bounds_x1, *bounds_x2),
            )
    plt.gca().invert_yaxis()
    if cbar:
        plt.colorbar()
    plt.axhline(0, color='k', linewidth=0.5)
    plt.axvline(0, color='k', linewidth=0.5)
    plt.xticks([-5, 0, +5])
    plt.yticks([0, 5, 10, 15, 20])
    return plt

def temp_diff(df: pd.DataFrame, stage: str, depth_idx: int) -> np.ndarray:
    after = df[f'{stage}_after_temp'].apply(
            lambda x: x[depth_idx]).values
    before = df[f'{stage}_before_temp'].apply(
            lambda x: x[depth_idx]).values
    return (after-before).copy()


