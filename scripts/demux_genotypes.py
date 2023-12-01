import gzip
import pysam
import click

from singlecell_snptools.io import read_genotypes


def geno_demux(bam_fn, cb_genotypes, output_fastq_pattern):
    genotypes = list(set(cb_genotypes.values()))
    output_fastq = {geno: gzip.open(output_fastq_pattern.format(geno=geno), 'wb') for geno in genotypes}
    with pysam.AlignmentFile(bam_fn) as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_secondary or aln.is_supplementary:
                continue
            else:
                cb = aln.get_tag('CB')
                umi = aln.get_tag('UB')
                if umi is None:
                    continue
                try:
                    geno = cb_genotypes[cb]
                except KeyError:
                    continue
                read_id = aln.query_name
                seq = aln.query_sequence
                qual = aln.qual
                output_fastq[geno].write(f'@{read_id}__{cb}__{umi}\n{seq}\n+\n{qual}\n'.encode())
    for fh in output_fastq:
        fh.close()



@click.command()
@click.option('-g', '--genotype-assignments-fn', required=True)
@click.option('-b', '--bam-fn', required=True)
@click.option('--pattern', required=True)
def main(genotype_assignments_fn, bam_fn, pattern):
    genotypes = read_genotypes(genotype_assignments_fn)
    geno_demux(bam_fn, genotypes, pattern)


if __name__ == '__main__':
    main()
