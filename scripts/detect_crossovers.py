import os

import numpy as np
import pandas as pd

from singlecell_snptools.io import get_co_markers_from_star_diploid_bam, read_genotypes, get_chrom_sizes_bam
from singlecell_snptools.crossovers import RigidHMM

from joblib import Parallel, delayed
import click


def create_training_data(co_markers, min_len, noise=1e-2, nexamples=1000):
    training_data = []
    for _, cb_co_markers in co_markers.groupby(['cb', 'chrom']):
        lr = cb_co_markers.lr.values
        if len(lr) >= min_len:
            if noise:
                lr = lr + np.random.normal(scale=noise, size=len(lr))
            training_data.append(lr.reshape(-1, 1))
        if len(training_data) >= nexamples:
            break
    return training_data


def train_rhmm(co_markers, rigidity=20, terminal_rigidity=8, training_nexamples=5000):
    X_fit = create_training_data(co_markers, min_len=terminal_rigidity, nexamples=training_nexamples)
    rhmm = RigidHMM(rigidity, term_rfactor=terminal_rigidity)
    rhmm.fit(X_fit, verbose=True)
    return rhmm


def detect_crossovers(co_markers, genotypes, chrom_sizes, rhmm, ref_name='col0'):
    invs = []
    for cb, cb_co_markers in co_markers.groupby('cb'):
        geno_states = [ref_name, genotypes[cb]]
        for chrom, m in cb_co_markers.groupby('chrom'):
            x = m.lr.values.reshape(-1, 1)
            starts, ends, start_marker, end_marker, states, logprobs, nmarkers = rhmm.predict(m.pos.values, x, chrom_sizes[chrom])
            for i, j, m, n, s, lp, nm in zip(starts, ends, start_marker, end_marker, states, logprobs, nmarkers):
                invs.append([chrom, i, j, m, n, cb, geno_states[s], lp, nm])
    invs = pd.DataFrame(invs, columns=['chrom', 'start', 'end', 'start_marker', 'end_marker', 'cb', 'geno', 'logprob', 'nmarkers'])
    invs = invs.sort_values(['chrom', 'start', 'cb'])
    return invs


def write_bed(invs, output_bed_fn):
    with open(output_bed_fn, 'w') as bed:
        for _, r in invs.iterrows():
            bed.write(f'{r.chrom:s}\t{r.start:d}\t{r.end:d}\t{r.cb}|{r.geno}|{r.nmarkers:d}|{r.start_marker}-{r.end_marker}\t{r.logprob}\t.\n')


@click.command()
@click.argument('bam_fns', nargs=-1)
@click.option('-g', '--genotype-assignments-fn', required=True)
@click.option('-m', '--model-fn', required=False, default=None)
@click.option('--output-model-fn', required=False, default=None)
@click.option('-o', '--output-bed-fn', required=True)
@click.option('--train-only/--train-predict', default=False)
@click.option('--rigidity', default=20)
@click.option('--terminal-rigidity', default=4)
@click.option('--training-ncells', default=1000)
@click.option('-p', '--processes', default=1)
def main(bam_fns, genotype_assignments_fn,
         model_fn, output_model_fn, output_bed_fn,
         train_only, rigidity, terminal_rigidity, training_ncells,
         processes):
    chrom_sizes = get_chrom_sizes_bam(bam_fns[0])
    genotypes = read_genotypes(genotype_assignments_fn)
    with Parallel(n_jobs=processes) as pool:
        co_markers = pool(
            delayed(get_co_markers_from_star_diploid_bam)(bam_fn)
            for bam_fn in bam_fns
        )
    co_markers = pd.concat(co_markers).sort_values(['cb', 'chrom', 'pos'])
    if model_fn is None:
        rhmm = train_rhmm(co_markers, rigidity, terminal_rigidity,
                          training_ncells * len(chrom_sizes))
        if output_model_fn is not None:
            rhmm.save(output_model_fn)
    else:
        rhmm = RigidHMM.load_model(model_fn)

    if not train_only:
        invs = detect_crossovers(co_markers, genotypes, chrom_sizes, rhmm)
        write_bed(invs, output_bed_fn)


if __name__ == '__main__':
    main()
