import matplotlib as mpl
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource

from skimage import measure

# Code attribution: https://stackoverflow.com/questions/56864378/how-to-light-and-shade-a-poly3dcollection
# Used via SE's default MIT license: https://meta.stackexchange.com/questions/272956/a-new-code-license-the-mit-this-time-with-attribution-required

rast=True

prefix = '_LOOCV'
lamb = 'AdaptInflate'
filt = 'Hur'
data_dir = 'Data/'
fig_dir = 'Figures_Depth/'
tps_est_dir = 'Estimates/'
tps_preds = pkl.load(open(f'{tps_est_dir}/TPS{prefix}_Preds_Block_{lamb}_{filt}.pkl',
        'rb'))
tps_masks = pkl.load(open(f'{tps_est_dir}/TPS{prefix}_Mask_Block_{lamb}_{filt}.pkl',
        'rb'))


tps_preds_mat = np.zeros((400, 100, 20))
mask_mat = np.zeros((400, 100, 20))
for depth in range(20):
    tps_preds_mat[:, :, 19-depth] = tps_preds[:, depth].reshape(400, 100)
    mask_mat[:, :, 19-depth] = tps_masks[:, depth].reshape(400, 100)

mask_mat = mask_mat.astype(bool)

fs = 16
df = 2
plt.rcParams['font.family'] = 'Liberation Serif'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['mathtext.it'] = 'serif:italic'
plt.rcParams['mathtext.bf'] = 'serif:bold'
plt.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams.update({'font.size': fs})



# --------------------------------------------------------------------------- #
# Restrict to first 10 days, less noisy
def plot_isosurface(init_elev, init_azim, output_fname): 
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=Axes3D.name)
    ax.view_init(init_elev, init_azim)


    lower_lim = -0.30
    upper_lim = +0.30

    tps_half = tps_preds_mat[:200, :, :]
    mask_half = mask_mat[:200, :, :]
    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    verts, faces, normals, values = measure.marching_cubes(tps_half, lower_lim,
	    mask=mask_half)
    mesh = Poly3DCollection(verts[faces], rasterized=rast)

    ls = LightSource(azdeg=225.0, altdeg=45.0)
    normalsarray = np.array([np.array((np.sum(normals[face[:], 0]/3),
	np.sum(normals[face[:], 1]/3), np.sum(normals[face[:],
	    2]/3))/np.sqrt(np.sum(normals[face[:], 0]/3)**2 +
		np.sum(normals[face[:], 1]/3)**2 + np.sum(normals[face[:],
		    2]/3)**2)) for face in faces])

    # Next this is more asthetic, but it prevents the shadows of the image being too dark. (linear interpolation to correct)
    min = np.min(ls.shade_normals(normalsarray, fraction=1.0)) # min shade value
    max = np.max(ls.shade_normals(normalsarray, fraction=1.0)) # max shade value
    diff = max-min
    newMin = 0.3
    newMax = 0.95
    newdiff = newMax-newMin

    # Using a constant color, put in desired RGB values here.
    colourRGB = np.array((0/255.0, 53.0/255.0, 107/255.0, 1.0))

    # The correct shading for shadows are now applied. Use the face normals and light orientation to generate a shading value and apply to the RGB colors for each face.
    rgbNew = np.array([colourRGB*(newMin + newdiff*((shade-min)/diff)) for shade in ls.shade_normals(normalsarray, fraction=1.0)])

    # Apply color to face
    mesh.set_facecolor(rgbNew)

    ax.add_collection3d(mesh)

    verts, faces, normals, values = measure.marching_cubes(tps_half, upper_lim,
	    mask=mask_half)
    mesh = Poly3DCollection(verts[faces], rasterized=rast)

    ls = LightSource(azdeg=225.0, altdeg=45.0)
    normalsarray = np.array([np.array((np.sum(normals[face[:], 0]/3),
	np.sum(normals[face[:], 1]/3), np.sum(normals[face[:],
	    2]/3))/np.sqrt(np.sum(normals[face[:], 0]/3)**2 +
		np.sum(normals[face[:], 1]/3)**2 + np.sum(normals[face[:],
		    2]/3)**2)) for face in faces])

    # Next this is more asthetic, but it prevents the shadows of the image being too dark. (linear interpolation to correct)
    min = np.min(ls.shade_normals(normalsarray, fraction=1.0)) # min shade value
    max = np.max(ls.shade_normals(normalsarray, fraction=1.0)) # max shade value
    diff = max-min
    newMin = 0.3
    newMax = 0.95
    newdiff = newMax-newMin

    # Using a constant color, put in desired RGB values here.
    colourRGB = np.array((255.0/255.0, 54.0/255.0, 57/255.0, 1.0))

    # The correct shading for shadows are now applied. Use the face normals and light orientation to generate a shading value and apply to the RGB colors for each face.
    rgbNew = np.array([colourRGB*(newMin + newdiff*((shade-min)/diff)) for shade in ls.shade_normals(normalsarray, fraction=1.0)])

    # Apply color to face
    mesh.set_facecolor(rgbNew)

    ax.add_collection3d(mesh)


    ax.set_xlabel(r"Days since TC passage ($\tau$)")
    ax.set_ylabel(r"Cross-track angle, degrees ($d$)")
    ax.set_zlabel(r"Pressure, dbars ($z$)")
    ax.set_xlim(0, 200)  
    ax.set_ylim(0, 100)  
    ax.set_zlim(0, 20)   

    ax.set_xticks([36, 91, 145, 199])
    ax.set_xticklabels([0, 3, 6, 9])
    ax.set_yticks([19, 50, 80])
    ax.set_yticklabels([-5, 0, 5])
    ax.set_zticks([0, 5, 10, 15, 19])
    ax.set_zticklabels([200, 150, 100, 50, 10])

    plt.savefig(output_fname,
            bbox_inches='tight',
            pad_inches=0,
            dpi=300,
            )
    return

init_elev1 = 21.0
init_azim1 = 144.1

# Second plot:
init_elev2 = 25.052
init_azim2 = 16.083

# Third plot:
init_elev3 = 18.557
init_azim3 = 176.412

plot_isosurface(
        init_elev1,
        init_azim1,
        'Figures_Isosurface/TPS_Isosurface_1.pdf',
        )

plot_isosurface(
        init_elev2,
        init_azim2,
        'Figures_Isosurface/TPS_Isosurface_2.pdf',
        )

plot_isosurface(
        init_elev3,
        init_azim3,
        'Figures_Isosurface/TPS_Isosurface_3.pdf',
        )

