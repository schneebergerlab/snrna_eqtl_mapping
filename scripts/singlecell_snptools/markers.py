import json
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import pysam
from joblib import Parallel, delayed

from .umi import umi_dedup_hap
from .io import ORGANELLAR_CONTIGS, get_chrom_sizes_bam, read_cb_whitelist

HOMOPOLYMERS = set(['AAAAAAAAAAAA', 'TTTTTTTTTTTT', 'CCCCCCCCCCCC', 'GGGGGGGGGGGG'])


def get_bin_umi_counts(bam, chrom, bin_idx, bin_size, whitelist_check):
    bin_umi_counts = defaultdict(lambda: defaultdict(Counter))
    bin_start = bin_idx * bin_size
    bin_end = bin_start + bin_size
    for aln in bam.fetch(chrom, bin_start, bin_end):
        if aln.is_secondary or aln.is_supplementary or aln.is_duplicate:
            continue
        # only consider fwd mapping read of properly paired 2xreads
        if aln.is_paired:
            if aln.is_reverse or not aln.is_proper_pair:
                continue
        if (aln.reference_start // bin_size) != bin_idx:
            continue
        cb = aln.get_tag('CB')
        if not whitelist_check(cb):
            continue
        hap = aln.get_tag('ha')
        umi = aln.get_tag('UB')
        if hap != 0:
            bin_umi_counts[cb][umi][hap] += 1

    bin_counts_dedup = defaultdict(Counter)
    for cb, umi_counts in bin_umi_counts.items():
        collapsed_umis = umi_dedup_hap(umi_counts)
        for hap_counts in collapsed_umis.values():
            if len(hap_counts) != 1:
                # UMI is ambiguous i.e. some reads map to hap1 and others to hap2
                continue
            else:
                hap = next(iter(hap_counts))
                bin_counts_dedup[cb][hap] += 1
    return bin_counts_dedup


def get_co_markers_from_star_diploid_bam(bam_fn, bin_size, whitelist=None, organellar_contigs=ORGANELLAR_CONTIGS):
    co_markers = defaultdict(dict)
    seen_cb = set()
    if whitelist is not None:
        whitelist_check = lambda cb: cb in whitelist
    else:
        whitelist_check = lambda cb: True
    with pysam.AlignmentFile(bam_fn) as bam:
        chrom_sizes = {k: bam.get_reference_length(k) for k in bam.references if k not in organellar_contigs}
        for chrom, cs in chrom_sizes.items():
            nbins = int(cs // bin_size + bool(cs % bin_size))
            chrom_co_markers = defaultdict(lambda: np.zeros(shape=(nbins, 2), dtype=np.uint32))
            for bin_idx in range(nbins):
                bin_co_markers = get_bin_umi_counts(
                    bam, chrom, bin_idx, bin_size, whitelist_check
                )
                for cb, hap_counts in bin_co_markers.items():
                    seen_cb.add(cb)
                    for hap, count in hap_counts.items():
                        chrom_co_markers[cb][bin_idx, hap - 1] = count
            for cb, m in chrom_co_markers.items():
                co_markers[cb][chrom] = m
    return co_markers, seen_cb


def estimate_bg_signal(co_markers, chrom_sizes, bin_size=25_000):
    bg_signal = {}
    for chrom, cs in chrom_sizes.items():
        nbins = int(cs // bin_size + bool(cs % bin_size))
        bg_signal[chrom] = np.zeros(shape=(nbins, 2), dtype=np.float64)

    for chrom_co_markers in co_markers.values():
        for chrom, m in chrom_co_markers.items():
            bg_signal[chrom] += m.astype(np.float64)
    for chrom, sig in bg_signal.items():
        bg_signal[chrom] = sig / sig.sum()
    return bg_signal


def estimate_and_subtract_background(co_markers, chrom_sizes, bin_size=25_000, window_size=1_000_000, max_bin_count=20):
    bg_signal = estimate_bg_signal(co_markers, chrom_sizes, bin_size)
    perc_contamination = {}
    co_markers_bg_subtracted = {}
    ws = window_size // bin_size
    for (fn_idx, cb), cb_co_markers in co_markers.items():
        nmarkers_tot = 0
        nmarkers_contam = 0
        for m in cb_co_markers.values():
            for i in range(len(m) - ws + 1):
                win = m[i: i + ws].sum(axis=0)
                nmarkers_tot += win.sum()
                nmarkers_contam += min(win)
        pc = nmarkers_contam / nmarkers_tot
        perc_contamination[(fn_idx, cb)] = pc
        bg_subtracted = {}
        for chrom, m in cb_co_markers.items():
            bg = pc * m.sum() * bg_signal[chrom]
            bg_subtracted[chrom] = np.minimum(
                np.round(np.maximum(m.astype(np.float64) - bg, 0)),
                max_bin_count
            ).astype(np.uint16)
        co_markers_bg_subtracted[(fn_idx, cb)] = bg_subtracted
    return perc_contamination, co_markers_bg_subtracted


def get_marker_ranges(co_markers, bin_size):
    marker_ranges = []
    for (fn_idx, cb), cb_co_markers in co_markers.items():
        for chrom, m in cb_co_markers.items():
            idx, = np.nonzero(m.sum(axis=1))
            left_pos, right_pos = idx[0] * bin_size, (idx[-1] + 1) * bin_size
            marker_ranges.append([fn_idx, cb, chrom, left_pos, right_pos])
    return pd.DataFrame(marker_ranges, columns=['fn_idx', 'cb', 'chrom', 'left_pos', 'right_pos'])


def co_markers_to_json(output_fn, co_markers, chrom_sizes, bin_size):
    co_markers_json_serialisable = {}
    marker_arr_sizes = {chrom: int(cs // bin_size + bool(cs % bin_size)) for chrom, cs, in chrom_sizes.items()}
    for (fn_idx, cb), cb_co_markers in co_markers.items():
        cb_id = f'{cb}-{fn_idx}'
        d = {}
        for chrom, arr in cb_co_markers.items():
            idx = np.nonzero(arr.ravel())[0]
            val = arr.ravel()[idx]
            d[chrom] = (idx.tolist(), val.tolist())
        co_markers_json_serialisable[cb_id] = d
    with open(output_fn, 'w') as o:
        return json.dump({
            'bin_size': bin_size,
            'shape': marker_arr_sizes,
            'data': co_markers_json_serialisable
        }, fp=o)


def load_co_markers_from_json(co_marker_json_fn, discard_fn_idx=False):
    with open(co_marker_json_fn) as f:
        co_marker_json = json.load(f)
    co_markers = {}
    bin_size = co_marker_json['bin_size']
    arr_shapes = co_marker_json['shape']
    for cb_id, cb_marker_idx in co_marker_json['data'].items():
        cb_co_markers = {}
        for chrom, (idx, val) in cb_marker_idx.items():
            m = np.zeros(shape=arr_shapes[chrom] * 2, dtype=np.uint16)
            m[idx] = val
            cb_co_markers[chrom] = m.reshape(-1, 2)
        cb, fn_idx = cb_id.split('-')
        co_markers[cb if discard_fn_idx else (int(fn_idx, cb))] = cb_co_markers
    return co_markers


def parallel_read_bams(bam_fns, bin_size, processes=1, whitelist_fn=None, whitelist=None, max_bin_count=20, bg_window_size=1_000_000):
    chrom_sizes = get_chrom_sizes_bam(bam_fns[0])
    if whitelist is None and whitelist_fn is not None:
        whitelist = read_cb_whitelist(whitelist_fn)
    with Parallel(n_jobs=min(processes, len(bam_fns)), backend='loky', verbose=True) as pool:
        res = pool(
            delayed(get_co_markers_from_star_diploid_bam)(bam_fn, bin_size, whitelist=whitelist)
            for bam_fn in bam_fns
        )
    co_markers = {}
    for f_idx, (f_co_markers, f_cb_whitelist) in enumerate(res):
        for cb in f_cb_whitelist:
            co_markers[(f_idx, cb)] = f_co_markers[cb]
    perc_contamination, co_markers = estimate_and_subtract_background(co_markers, chrom_sizes, bin_size, bg_window_size, max_bin_count)
    marker_ranges = get_marker_ranges(co_markers, bin_size)
    return co_markers, marker_ranges, perc_contamination, chrom_sizes
