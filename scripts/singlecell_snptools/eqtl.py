import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import convolve
import pandas as pd

from .io import parse_gene_positions

LOD_THRESHOLD = 3.0
REL_LOD_THRESHOLD = 0.1
FDR_THRESHOLD = 0.05
REL_PROMININENCE = 0.25
LOD_DROP = 1.5
MIN_DIST = 3_000_000
BIN_SIZE = 25_000


def find_lod_peaks(eqtl_data, lod_threshold=LOD_THRESHOLD, rel_lod_threshold=REL_LOD_THRESHOLD,
                   fdr_threshold=FDR_THRESHOLD,
                   rel_prominence=REL_PROMININENCE, lod_drop=LOD_DROP,
                   min_dist=MIN_DIST, bin_size=BIN_SIZE,
                   cell_type_specific=False, fdr_col=None):
    min_dist = min_dist // bin_size
    peaks = []
    for _, chrom_eqtl_data in eqtl_data.groupby('chrom'):
        chrom_eqtl_data = chrom_eqtl_data.sort_values('pos')
        lod = chrom_eqtl_data.lod_score.values if not cell_type_specific else np.maximum(chrom_eqtl_data.cluster_lod_score.values, chrom_eqtl_data.lod_score.values)
        lod[lod < lod_threshold] = -np.inf
        if fdr_col is None:
            fdr = chrom_eqtl_data.lrt_pval.values if not cell_type_specific else np.minimum(chrom_eqtl_data.cluster_lrt_pval.values, chrom_eqtl_data.lrt_pval.values) * 2
        else:
            fdr = chrom_eqtl_data[fdr_col].values
        peak_idx, _ = find_peaks(
            np.insert(lod, [0, len(lod)], [-np.inf, -np.inf]), # insert infs at beginning and end of chromosome to prevent prominence edge effects
            height=max(lod_threshold, max(lod) * rel_lod_threshold),
            prominence=max(lod) * rel_prominence,
            distance=min_dist
        )
        peak_idx = peak_idx - 1
        peak_idx = [i for i in peak_idx if fdr[i] < fdr_threshold]
        l_ci_idx = 0
        r_ci_idx = len(chrom_eqtl_data) - 1
        for p_idx in peak_idx:
            drop_thresh = lod[p_idx] - lod_drop
            for r_ci_idx in range(p_idx, len(lod)):
                if lod[r_ci_idx] <= drop_thresh:
                    break
            for l_ci_idx in reversed(range(0, p_idx)):
                if lod[l_ci_idx] <= drop_thresh:
                    break
            eqtl_peak = chrom_eqtl_data.iloc[p_idx].copy()
            eqtl_peak['start_eqtl'] = chrom_eqtl_data.iloc[l_ci_idx].pos
            eqtl_peak['end_eqtl'] = chrom_eqtl_data.iloc[r_ci_idx].pos + bin_size
            peaks.append(eqtl_peak)
    peaks = pd.DataFrame(
        peaks,
        columns=eqtl_data.columns.tolist() + ['start_eqtl', 'end_eqtl'],
    )
    peaks['pos_eqtl'] = peaks.pop('pos').astype(int)
    peaks['start_eqtl'] = peaks.start_eqtl.astype(int)
    peaks['end_eqtl'] = peaks.end_eqtl.astype(int)
    return peaks


def call_lod_peaks(eqtl_res, gtf_fn, cis_range=2e6, **kwargs):
    sig_eqtl = eqtl_res.groupby('gene_id', as_index=False).apply(find_lod_peaks, **kwargs).reset_index(drop=True)
    gene_locs = parse_gene_positions(gtf_fn)
    gene_locs = pd.DataFrame.from_dict(gene_locs, orient='index', columns=['chrom', 'pos_gene'])
    sig_eqtl = pd.merge(
        sig_eqtl,
        gene_locs,
        left_on='gene_id',
        right_index=True,
        how='left',
        suffixes=['_eqtl', '_gene']
    )
    sig_eqtl = sig_eqtl.dropna(subset=['chrom_gene'])
    sig_eqtl['is_cis'] = sig_eqtl.eval('chrom_eqtl == chrom_gene & (start_eqtl - @cis_range) <= pos_gene <= (end_eqtl + @cis_range)')
    return sig_eqtl