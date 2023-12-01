import numpy as np
from scipy.signal import find_peaks
import pandas as pd


def find_lod_peaks(eqtl_data, nlog10fdr_threshold=3, lod_drop=1.5, min_dist=250):
    peaks = []
    for _, chrom_eqtl_data in eqtl_data.groupby('chrom'):
        chrom_eqtl_data = chrom_eqtl_data.sort_values('start')
        lod = chrom_eqtl_data.lod_score.values
        fdr = -np.log10(chrom_eqtl_data.lrt_fdr.values)
        peak_idx, peak_stats = find_peaks(
            fdr,
            height=nlog10fdr_threshold,
            prominence=nlog10fdr_threshold / 2,
            distance=min_dist
        )
        for p_idx, h, lb, rb in zip(peak_idx, peak_stats['peak_heights'], peak_stats['left_bases'], peak_stats['right_bases']):
            drop_thresh = h - lod_drop
            for r_ci_idx in range(p_idx, rb):
                if lod[r_ci_idx] <= drop_thresh:
                    break
            for l_ci_idx in reversed(range(lb, p_idx)):
                if lod[l_ci_idx] <= drop_thresh:
                    break
            eqtl_peak = chrom_eqtl_data.iloc[p_idx].copy()
            eqtl_peak['pos_eqtl'] = (eqtl_peak.start + eqtl_peak.end) // 2
            eqtl_peak['start'] = chrom_eqtl_data.iloc[l_ci_idx].start
            eqtl_peak['end'] = chrom_eqtl_data.iloc[r_ci_idx].end
            peaks.append(eqtl_peak)
    peaks =  pd.DataFrame(
        peaks,
        columns=eqtl_data.columns.tolist() + ['pos_eqtl'],
    )
    peaks['pos_eqtl'] = peaks.pos_eqtl.astype(int)
    peaks['start'] = peaks.start.astype(int)
    peaks['end'] = peaks.end.astype(int)
    peaks = peaks.loc[:, ['chrom', 'pos_eqtl'] + peaks.columns.tolist()[1:-1]]
    return peaks