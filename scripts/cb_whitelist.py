import os
from collections import Counter
import pandas as pd
import click


def read_index(index_fn):
    with open(index_fn) as f:
        idx = [r.strip() for r in f.readlines()]
    return idx


def quick_count(mtx_fn, feature_fn, barcode_fn, mt_id='ATMG', cp_id='ATCG'):

    barcodes = read_index(barcode_fn)
    features = read_index(feature_fn)
    feature_is_mt = set([i for i, gene_id in enumerate(features) if gene_id.startswith(mt_id)])
    feature_is_cp = set([i for i, gene_id in enumerate(features) if gene_id.startswith(cp_id)])
    
    cb_umi_count = Counter()
    cb_mt_count = Counter()
    cb_cp_count = Counter()

    with open(mtx_fn) as mtx:
        while True:
            line = next(mtx)
            if line.startswith('%'):
                continue
            else:
                nrow, ncol, nonzero = line.split()
                break
        for line in mtx:
            row, col, count = line.split()
            row = int(row) - 1
            col = int(col) - 1
            count = int(count)
            cb = barcodes[col]
            cb_umi_count[cb] += count
            if row in feature_is_mt:
                cb_mt_count[cb] += count
            elif row in feature_is_cp:
                cb_cp_count[cb] += count

    cb_stats = []
    for cb in barcodes:
        tot = cb_umi_count[cb]
        if tot:
            cb_stats.append([cb, tot, cb_mt_count[cb] / tot * 100, cb_cp_count[cb] / tot * 100])
        else:
            cb_stats.append([cb, 0, 0.0, 0.0])
    return pd.DataFrame(cb_stats, columns=['cb', 'umis', 'perc_mt', 'perc_cp'])


@click.command()
@click.option('-s', '--starsolo-dir', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('--min-umis-per-barcode', default=100)
@click.option('--max-perc-mito', default=5)
@click.option('--max-perc-chloro', default=5)
def main(starsolo_dir, output_fn, min_umis_per_barcode, max_perc_mito, max_perc_chloro):
    cb_stats = quick_count(
        os.path.join(starsolo_dir, 'Gene/raw/matrix.mtx'),
        os.path.join(starsolo_dir, 'Gene/raw/features.tsv'),
        os.path.join(starsolo_dir, 'Gene/raw/barcodes.tsv'),
    )
    cb_stats = cb_stats.query('umis >= @min_umis_per_barcode & perc_mt < @max_perc_mito & perc_cp < @max_perc_chloro')
    with open(output_fn, 'w') as f:
        for cb in cb_stats.cb.values:
            f.write(f'{cb}\n')
            
if __name__ == '__main__':
    main()