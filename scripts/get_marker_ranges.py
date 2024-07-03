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


@click.command()
@click.argument('bam_fns', nargs=-1)
@click.option('-o', '--output-prefix', required=False)
@click.option('-w', '--whitelist-fn', required=False, default=None)
@click.option('--bin-size', default=25_000)
@click.option('-p', '--processes', default=1)
def main(bam_fns, output_prefix, whitelist_fn, bin_size, processes):
    print(f'{len(bam_fns)} bams provided')
    bam_fns = sorted(bam_fns)
    print('loading crossover markers')
    co_markers, marker_ranges, perc_contamination, chrom_sizes = parallel_read_bams(
        bam_fns,
        bin_size,
        processes=processes,
        whitelist_fn=whitelist_fn,
    )
    co_markers_to_json(f'{output_prefix}.co_markers.json', co_markers, chrom_sizes, bin_size)
    marker_ranges.to_csv(f'{output_prefix}.marker_ranges.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
