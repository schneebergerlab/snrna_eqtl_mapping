import os
import re
from copy import copy
from collections import defaultdict, Counter
import itertools as it
from scipy import sparse
from scipy.io import mmwrite
import numpy as np

import pysam
from joblib import Parallel, delayed
import click

from singlecell_snptools.io import read_cb_whitelist, parse_gtf
from singlecell_snptools.umi import umi_dedup
from singlecell_snptools.markers import HOMOPOLYMERS


def strand_filter(filt_type, gene_strand):
    gene_is_rev = gene_strand == '-'
    if filt_type == 'fwd':
        def _strand_filter(aln):
            return aln.is_reverse == gene_is_rev
    elif filt_type == 'rev':
        def _strand_filter(aln):
            return aln.is_reverse != gene_is_rev
    else:
        def _strand_filter(aln):
            return True
    return _strand_filter


def generate_scrna_count_matrix(genes, cb_whitelist, bam_fn, min_mapq=255, strand_filt_type='fwd'):
    gene_counts = sparse.dok_matrix((len(genes), len(cb_whitelist)), dtype=np.int16)
    cb_idx = {cb: i for i, cb in enumerate(cb_whitelist)}
    with pysam.AlignmentFile(bam_fn) as bam:
        for gene_idx, (chrom, gene_strand, exons) in enumerate(genes):
            sf = strand_filter(strand_filt_type, gene_strand)
            umi_counts = defaultdict(Counter)
            for start, end in exons:
                for aln in bam.fetch(chrom, start, end):
                    if aln.mapq >= min_mapq and sf(aln):
                        cb = aln.get_tag('CB')
                        umi = aln.get_tag('UR')
                        if umi.count('N') or umi in HOMOPOLYMERS:
                            continue
                        else:
                            try:
                                c_i = cb_idx[cb]
                            except KeyError:
                                continue
                            umi_counts[c_i][umi] += 1

            for c_i, cb_umi_counts in umi_counts.items():
                cb_umi_counts = umi_dedup(cb_umi_counts)
                gene_counts[gene_idx, c_i] = len(cb_umi_counts)

    return gene_counts.tocsc()


def chunk_data(data, chunk_size=100):
    d_it = iter(data)
    for i in range(0, len(data), chunk_size):
        yield list(it.islice(d_it, chunk_size))


def parallel_scrna_gene_counter(gtf_fn, cb_whitelist_fn, bam_fn, min_mapq=255, strand='fwd', quant_type='gene', njobs=1):
    cb_whitelist = read_cb_whitelist(cb_whitelist_fn)
    gene_ids, gene_locs = parse_gtf(gtf_fn, quant_type=quant_type)
    with Parallel(njobs, verbose=100) as pool:
        matrices = pool(
            delayed(generate_scrna_count_matrix)(
                gene_chunk, cb_whitelist, bam_fn, min_mapq, strand_filt_type=strand
            ) for gene_chunk in list(chunk_data(gene_locs))
        )
    return gene_ids, cb_whitelist, sparse.vstack(matrices)


def write_features(feats_list, output_fn):
    with open(output_fn, 'w') as f:
        for feat in feats_list:
            f.write(f'{feat}\n')


@click.command()
@click.option('-o', '--output-dir', required=True)
@click.option('-g', '--gtf-fn', required=True)
@click.option('-w', '--cb-whitelist-fn', required=True)
@click.option('-b', '--bam-fn', required=True)
@click.option('-q', '--min-mapq', required=False, default=255)
@click.option('-s', '--strand', type=click.Choice(['fwd', 'rev', 'both']), default='fwd')
@click.option('--quant-type', type=click.Choice(['gene', 'gene_full']), default='gene')
@click.option('-n', '--processes', required=False, default=1)
def main(output_dir, gtf_fn, cb_whitelist_fn, bam_fn, min_mapq, strand, quant_type, processes):
    gene_ids, cb_whitelist, mtx = parallel_scrna_gene_counter(
        gtf_fn, cb_whitelist_fn, bam_fn, min_mapq, strand, quant_type, njobs=processes
    )
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    cb_fn = os.path.join(output_dir, 'barcodes.tsv')
    feat_fn = os.path.join(output_dir, 'features.tsv')
    mtx_fn = os.path.join(output_dir, 'matrix.mtx')

    write_features(cb_whitelist, cb_fn)
    write_features(gene_ids, feat_fn)
    mmwrite(mtx_fn, mtx)


if __name__ == '__main__':
    main()