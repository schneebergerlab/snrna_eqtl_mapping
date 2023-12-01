import numpy as np
import pandas as pd
from scipy.io import mmread


def read_index(idx_fn):
    with open(idx_fn) as f:
        idx = [r.strip().split('\t')[0] for r in f.readlines()]
    return np.array(idx)


def lognormalise_exprs(mtx):
    mtx = (mtx / mtx.sum(0)) * 10_000
    lm_mtx = np.log2(mtx)
    lm_mtx[mtx < 1] = 0
    return lm_mtx


def read_mmformat(mtx_fn, barcode_fn, feature_fn, bc_filt=200, feat_filt=3):
    mm = mmread(mtx_fn).tocsr()
    bc_mask = np.array(mm.sum(0) > bc_filt).ravel()
    mm = mm[:, bc_mask]
    feat_mask = np.array(mm.sum(1) > feat_filt).ravel()
    mm = mm[feat_mask]
    mm = mm.todense().astype(np.float32)
    ln_mm = lognormalise_exprs(mm)
    barcodes = read_index(barcode_fn)[bc_mask]
    features = read_index(feature_fn)[feat_mask]
    mm = pd.DataFrame(mm, columns=barcodes, index=features)
    ln_mm = pd.DataFrame(ln_mm, columns=barcodes, index=features)
    return ln_mm


def get_coverage(mtx_fn, barcode_fn):
    mm = mmread(mtx_fn).tocsr()
    cov = np.array(mm.sum(0)).ravel()
    return pd.Series(np.log2(cov), index=read_index(barcode_fn), name='coverage')


def get_n_expressed(mtx_fn, barcode_fn):
    mm = mmread(mtx_fn).tocsr()
    n_exprs = np.array((mm > 0).sum(0)).ravel()
    return pd.Series(np.log2(n_exprs), index=read_index(barcode_fn), name='expressed_genes')