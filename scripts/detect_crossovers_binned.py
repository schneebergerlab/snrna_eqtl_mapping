import os
import json
from copy import copy
import itertools as it
from collections import defaultdict, Counter

os.environ['OPENBLAS_NUM_THREADS'] = '1'

import numpy as np
import pandas as pd
import pomegranate as pm
import pysam

from joblib import Parallel, delayed
import click

from singlecell_snptools.markers import parallel_read_bams, co_markers_to_json
from singlecell_snptools.rhmm import create_training_data, RigidHMM


def train_rhmm(co_markers, bin_size=25_000, cm_per_mb=4.5,
               min_seg_length=2_500_000, min_terminal_seg_length=100_000,
               training_nexamples=100):
    rigidity = int(min_seg_length // bin_size)
    terminal_rigidity = min_terminal_seg_length = int(min_terminal_seg_length // bin_size)
    co_prob = cm_per_mb * (bin_size / 1e8)
    X_fit = create_training_data(co_markers, nexamples=training_nexamples)
    print(f'Training with rigidity {rigidity} and terminal rigidity {terminal_rigidity}')
    rhmm = RigidHMM(rigidity, term_rfactor=terminal_rigidity, trans_prob=co_prob)
    rhmm.fit(X_fit, verbose=True, stop_threshold=1, multiple_check_input=False)
    return rhmm


def detect_crossovers(co_markers, chrom_sizes, rhmm, bin_size):
    co_probs = {}
    cb_stats = []
    expected_nbins = {c: int(cs // bin_size + bool(cs % bin_size)) for c, cs in chrom_sizes.items()}
    for (f_idx, cb), cb_co_markers in co_markers.items():
        cb_probs = {}
        cb_classification_score = []
        cb_nmarkers = 0
        for chrom, x in cb_co_markers.items():
            assert len(x) == expected_nbins[chrom], f'{cb} {chrom} {len(x)} {expected_nbins[chrom]}'
            x = x.astype(np.float32)
            cb_nmarkers += x.sum()
            states, probs = rhmm.predict(x)
            cls_score = (x[:, 0] * (1 - probs)).sum() + (x[:, 1] * probs).sum()
            cb_classification_score.append(cls_score)
            assert len(probs) == expected_nbins[chrom], f'{cb} {chrom} {len(x)} {len(probs)} {expected_nbins[chrom]}'
            cb_probs[chrom] = probs

        cb_classification_score = sum(cb_classification_score) / cb_nmarkers
        cb_stats.append([f_idx, cb, cb_nmarkers, cb_classification_score])

        cb_probs_concat = []
        for chrom in chrom_sizes:
            try:
                cb_probs_concat.append(cb_probs[chrom])
            except KeyError:
                cb_probs_concat.append(np.full(shape=expected_nbins[chrom], fill_value=0.5))
        co_probs[(f_idx, cb)] = np.concatenate(cb_probs_concat)

    cb_stats = pd.DataFrame(cb_stats, columns=['fn_idx', 'cb', 'nmarkers', 'classification_score'])
    co_probs = pd.DataFrame.from_dict(co_probs, orient='index')
    co_probs.columns = pd.MultiIndex.from_tuples(zip(
        np.repeat(list(expected_nbins.keys()), list(expected_nbins.values())),
        np.concatenate([np.arange(0, e) for e in expected_nbins.values()]) * bin_size
    ), names=['chrom', 'pos'])
    co_probs.index = pd.MultiIndex.from_tuples(co_probs.index.values, names=['fn_idx', 'cb'])
    return cb_stats, co_probs


def chunk_data(data, chunk_size):
    d_it = iter(data)
    for i in range(0, len(data), chunk_size):
        yield {k: data[k] for k in it.islice(d_it, chunk_size)}


def parallel_predict_crossovers(co_markers, chrom_sizes, rhmm, bin_size, processes=1):
    chunk_size = min(len(co_markers) // processes, 100)
    with Parallel(n_jobs=processes, backend='loky', verbose=True) as pool:
        res = pool(
            delayed(detect_crossovers)(co_chunk, chrom_sizes, rhmm, bin_size)
            for co_chunk in list(chunk_data(co_markers, chunk_size))
        )
    cb_stats, co_probs = zip(*res)
    cb_stats = pd.concat(cb_stats, axis=0)
    co_probs = pd.concat(co_probs, axis=0)
    return cb_stats, co_probs


@click.command()
@click.argument('bam_fns', nargs=-1)
@click.option('-m', '--model-fn', required=False, default=None)
@click.option('--output-model-fn', required=False, default=None)
@click.option('--output-marker-json-fn', required=False, default=None)
@click.option('-w', '--whitelist-fn', required=False, default=None)
@click.option('-o', '--output-prefix', required=False)
@click.option('--train-only/--train-predict', default=False)
@click.option('--bin-size', default=25_000)
@click.option('--cm-per-mb', default=4.5)
@click.option('--min-segment-size', default=2_500_000)
@click.option('--min-terminal-segment-size', default=100_000)
@click.option('--training-ncells', default=500)
@click.option('-p', '--processes', default=1)
def main(bam_fns, model_fn, output_model_fn, output_marker_json_fn, whitelist_fn, output_prefix,
         train_only, bin_size, cm_per_mb, min_segment_size, min_terminal_segment_size,
         training_ncells, processes):
    print(f'{len(bam_fns)} bams provided')
    bam_fns = sorted(bam_fns)
    training_ncells = training_ncells // len(bam_fns)
    if model_fn is None and train_only:
        print(f'selecting {training_ncells} cells from each bam')
    print('loading crossover markers')
    co_markers, marker_ranges, perc_contamination, chrom_sizes = parallel_read_bams(
        bam_fns,
        bin_size,
        processes=processes,
        whitelist_fn=whitelist_fn,
    )

    if output_marker_json_fn is not None:
        co_markers_to_json(output_marker_json_fn, co_markers, chrom_sizes, bin_size)
    
    if model_fn is None:
        rhmm = train_rhmm(co_markers,
                          bin_size=bin_size,
                          cm_per_mb=cm_per_mb,
                          min_seg_length=min_segment_size,
                          min_terminal_seg_length=min_terminal_segment_size,
                          training_nexamples=training_ncells * len(chrom_sizes))
        if output_model_fn is not None:
            rhmm.save(output_model_fn)
    else:
        rhmm = RigidHMM.load_model(model_fn)

    if not train_only:
        print(f'predicting crossovers for {len(co_markers)} cell barcodes')
        cb_stats, probs = parallel_predict_crossovers(co_markers, chrom_sizes, rhmm, bin_size=bin_size, processes=processes)
        cb_stats['fn'] = cb_stats.fn_idx.map(dict(enumerate(bam_fns)))
        cb_stats.to_csv(f'{output_prefix}.cb_stats.tsv', sep='\t', index=False)
        marker_ranges.to_csv(f'{output_prefix}.marker_ranges.tsv', sep='\t', index=False)
        probs.to_csv(f'{output_prefix}.co_probs.tsv', sep='\t')


if __name__ == '__main__':
    main()
