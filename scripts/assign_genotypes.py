import os
import random
from collections import defaultdict, Counter

os.environ['OPENBLAS_NUM_THREADS'] = '1'

import numpy as np
from scipy.io import mmread
import pandas as pd

import pysam
from joblib import Parallel, delayed
import click

from singlecell_snptools.io import parse_snps, parse_cellsnplite


def iter_cb(snp_counts):
    for cb in set(snp_counts.columns.get_level_values(0)):
        cb_snp_counts = snp_counts[cb]
        cb_snp_counts = cb_snp_counts.iloc[
            np.unique(np.concatenate([
                cb_snp_counts.ref_count.values.sp_index.indices,
                cb_snp_counts.alt_count.values.sp_index.indices,
            ]))
        ]
        cb_snp_counts = cb_snp_counts.sparse.to_dense().reset_index()
        if len(cb_snp_counts):
            yield cb, cb_snp_counts


def update_probs(probs, sample_markers):
    marker_agg = Counter()
    for _, m in sample_markers.iterrows():
        n_ref = len(m.ref_accs)
        n_alt = len(m.alt_accs)
        for acc in m.ref_accs:
            marker_agg[acc] += (probs[acc] * (bool(m.ref_count))) / n_ref
        for acc in m.alt_accs:
            marker_agg[acc] += (probs[acc] * (bool(m.alt_count))) / n_alt

    agg_sum = sum(marker_agg.values())
    return {acc: marker_agg[acc] / agg_sum for acc in probs}


def assign_genotype_with_em(sample_markers, accessions, max_iter=1000, min_delta=1e-2):
    n_accs = len(accessions)
    probs = {acc: 1 / n_accs for acc in accessions}
    for _ in range(max_iter):
        prev_probs = probs
        probs = update_probs(probs, sample_markers)
        delta = sum(abs(prev_probs[g] - probs[g]) for g in accessions)
        if delta < min_delta:
            break
    return probs


def bootstrap_genotype_assignment(cb, cb_snp_counts, accessions,
                                  em_max_iter, em_min_delta,
                                  bootstrap_sample_size, n_bootstraps):
    prob_boots = []
    cb_snp_counts = cb_snp_counts.query('alt_count > ref_count')
    if not len(cb_snp_counts):
        return cb, None, 1, 0.0
    bootstrap_sample_size = min(len(cb_snp_counts), bootstrap_sample_size)
    for _ in range(n_bootstraps):
        samp = cb_snp_counts.sample(bootstrap_sample_size, replace=True)
        p = assign_genotype_with_em(samp, accessions)
        prob_boots.append(p)
    prob_boots = pd.DataFrame.from_dict(prob_boots)
    prob_means = prob_boots.mean()
    assignment = prob_means.idxmax()
    pm = prob_means.loc[assignment]
    return cb, assignment, len(cb_snp_counts), -np.log10(1 - pm)


def parallel_assign_genotypes(snp_counts, accessions,
                              em_max_iter, em_min_delta,
                              bootstrap_sample_size, n_bootstraps,
                              processes):
    if len(accessions) > 1:
        n_cb = len(snp_counts.columns)
        kwargs = dict(
            em_max_iter=em_max_iter,
            em_min_delta=em_min_delta,
            bootstrap_sample_size=bootstrap_sample_size,
            n_bootstraps=n_bootstraps
        )
        assignments = Parallel(processes, verbose=11)(
            delayed(bootstrap_genotype_assignment)(cb, cb_snp_counts, accessions, **kwargs)
            for cb, cb_snp_counts in iter_cb(snp_counts)
        )
    else:
        # don't need to do anything! dummy assignment
        assignments = [
            [cb, accessions[0], len(cb_snp_counts), 1000.]
            for cb, cb_snp_counts in iter_cb(snp_counts)
        ]
    return pd.DataFrame(
        assignments,
        columns=['cell_barcode', 'assignment', 'nmarkers', 'logprob']
    ).query('assignment == assignment') # removes NoneType assignments


@click.command()
@click.option('-c', '--cellsnplite-dir', required=True)
@click.option('-v', '--vcf-fn', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('--em-max-iter', default=1000, type=int)
@click.option('--em-min-delta', default=1e-2, type=float)
@click.option('--bootstrap-sample-size', default=1000, type=int)
@click.option('--n-bootstraps', default=25, type=int)
@click.option('-p', '--processes', default=1, type=int)
def main(cellsnplite_dir, vcf_fn, output_fn,
         em_max_iter, em_min_delta,
         bootstrap_sample_size, n_bootstraps,
         processes):
    barcode_fn = os.path.join(cellsnplite_dir, 'cellSNP.samples.tsv')
    ad_fn = os.path.join(cellsnplite_dir, 'cellSNP.tag.AD.mtx')
    dp_fn = os.path.join(cellsnplite_dir, 'cellSNP.tag.DP.mtx')
    snp_counts, accessions = parse_cellsnplite(vcf_fn, barcode_fn, ad_fn, dp_fn)
    genotype_assignments = parallel_assign_genotypes(
        snp_counts, accessions,
        em_max_iter, em_min_delta,
        bootstrap_sample_size, n_bootstraps,
        processes
    )
    genotype_assignments.to_csv(output_fn, sep='\t', index=False)

    
if __name__ == '__main__':
    main()