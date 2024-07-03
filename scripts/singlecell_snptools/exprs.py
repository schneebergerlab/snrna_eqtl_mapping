import os
import numpy as np
import pandas as pd
from scipy.io import mmread


def read_index(idx_fn):
    with open(idx_fn) as f:
        idx = [r.strip().split('\t')[0] for r in f.readlines()]
    return np.array(idx)


def lognormalise_exprs(mtx, alpha=0.25):
    umi_counts = mtx.sum(0)
    sf = (umi_counts / np.mean(umi_counts))
    ln_mtx = np.log2(mtx / sf + (1.0 / (4 * alpha)))
    sf = sf.rename('size_factors')
    return ln_mtx, sf


def read_mmformat(mtx_dirs, cb_whitelist, feat_filt=10,
                  mtx_fn='matrix.mtx',
                  barcode_fn='barcodes.tsv',
                  feature_fn='features.tsv'):
    mtxs = []
    for mtx_dir in mtx_dirs:
        mtx = mmread(os.path.join(mtx_dir, mtx_fn)).tocsr()
        bc_mask = np.array(mtx.sum(axis=0) > 0).ravel()
        mtx = mtx[:, bc_mask]
        barcodes = read_index(os.path.join(mtx_dir, barcode_fn))[bc_mask]
        features = read_index(os.path.join(mtx_dir, feature_fn))
        mtx = pd.DataFrame.sparse.from_spmatrix(mtx, columns=barcodes, index=features)
        mtx = mtx.loc[:, mtx.columns.isin(cb_whitelist)]
        mtxs.append(mtx)
    mtx = pd.concat(mtxs, axis=1).fillna(0)
    mtx = mtx.loc[np.count_nonzero(mtx, axis=1) > feat_filt]
    mtx = mtx.sparse.to_dense().astype(np.float32)
    return mtx
