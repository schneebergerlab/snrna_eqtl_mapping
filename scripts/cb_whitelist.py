import os
import click
from singlecell_snptools.exprs import read_mmformat


@click.command()
@click.option('-s', '--starsolo-dir', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('--minimum-reads-per-barcode', default=200)
def main(starsolo_dir, output_fn, minimum_reads_per_barcode):
    exprs = read_mmformat(
        os.path.join(starsolo_dir, 'Gene/raw/matrix.mtx'),
        os.path.join(starsolo_dir, 'Gene/raw/barcodes.tsv'),
        os.path.join(starsolo_dir, 'Gene/raw/features.tsv'),
        bc_filt=minimum_reads_per_barcode
    )
    with open(output_fn, 'w') as f:
        for cb in exprs.columns:
            f.write(f'{cb}\n')
            
if __name__ == '__main__':
    main()