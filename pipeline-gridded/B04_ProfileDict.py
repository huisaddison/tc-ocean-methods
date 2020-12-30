import h5py
from tools import create_ArgoProfileID
import scipy.io
import pickle as pkl
import argparse

parser = argparse.ArgumentParser(description='Create profile dictionaries')
parser.add_argument('--integrated', dest='mode', action='store_const',
                    const='integrated', default='gridded',
                    help='process assuming integrated heat content (default: gridded)')
args = parser.parse_args()

for fi, val in [('ProfFiltered', 'Prof'), ('Res', 'Res')]:
    mat = h5py.File(f'./Data/'
        f'gridTemp{fi}_8_20.mat')

    vals = mat[f'gridTemp{val}'].value
    if args.mode == 'gridded':
        vals = vals.T
    dict_ = {
        create_ArgoProfileID(arr_fid[0], arr_cn[0]) : arr_rid
        for arr_fid, arr_cn, arr_rid
        in zip(mat['profFloatIDAggrSel'].value,
               mat['profCycleNumberAggrSel'].value,
               vals,
              )
        }

    pkl.dump(dict_, open(f'./Data/Hurricane{val}Dict.pkl', 'wb'))

