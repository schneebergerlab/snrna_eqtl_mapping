import re
from collections import defaultdict
from statistics import median_low

import numpy as np
import pandas as pd
from scipy.io import mmread

import pysam

ORGANELLAR_CONTIGS = set(['ChrM', 'ChrC'])


def get_chrom_sizes_bam(bam_fn, organellar_contigs=ORGANELLAR_CONTIGS):
    with pysam.AlignmentFile(bam_fn) as bam:
        chrom_sizes = {k: bam.get_reference_length(k) for k in bam.references if k not in organellar_contigs}
    return chrom_sizes


def parse_snps(vcf_fn):
    with pysam.VariantFile(vcf_fn) as vcf:
        accessions = list(vcf.header.samples)
        parsed_snps = []
        for snp in vcf.fetch():
            acc_alleles = defaultdict(list)
            acc_alleles[snp.ref] = []
            for acc, var in snp.samples.items():
                alt = var.alleles[0]
                acc_alleles[alt].append(acc)
            try:
                alt = [a for a in acc_alleles if a != snp.ref][0]
            except IndexError:
                alt = None
            parsed_snps.append([snp.contig, snp.pos, snp.ref, alt, tuple(acc_alleles[snp.ref]), tuple(acc_alleles[alt])])
        return parsed_snps, accessions


def read_cb_whitelist(barcode_fn):
    with open(barcode_fn) as f:
        cb_whitelist = [cb.strip().split('\t')[0] for cb in f.readlines()]
    return cb_whitelist


def read_genotypes(genotypes_fn):
    genotypes = pd.read_csv(genotypes_fn, sep='\t', dtype={'assignment': 'category'})
    return genotypes.set_index('cell_barcode').assignment.to_dict()


def parse_cellsnplite(vcf_fn, barcode_fn, ad_fn, dp_fn):
    alt_mm = mmread(ad_fn)
    dp_mm = mmread(dp_fn)
    ref_mm = dp_mm - alt_mm
    del dp_mm
    barcodes = read_cb_whitelist(barcode_fn)
    snps, accessions = parse_snps(vcf_fn)
    snps = pd.MultiIndex.from_tuples(snps, names=['chrom', 'pos', 'ref', 'alt', 'ref_accs', 'alt_accs'])
    alt_mm = pd.DataFrame.sparse.from_spmatrix(alt_mm, columns=barcodes, index=snps)
    ref_mm = pd.DataFrame.sparse.from_spmatrix(ref_mm, columns=barcodes, index=snps)
    counts = pd.concat({'ref_count': ref_mm, 'alt_count': alt_mm}, axis=1)
    counts.columns = counts.columns.reorder_levels([1, 0])
    return counts, accessions


def get_gtf_attribute(gtf_record, attribute):
    try:
        attr = re.search(f'{attribute} "(.+?)";', gtf_record[8]).group(1)
    except AttributeError:
        raise ValueError(
            f'Could not parse attribute {attribute} '
            f'from GTF with feature type {record[2]}'
        )
    return attr


def parse_gene_positions(gtf_fn):
    gtf_records = {}
    with open(gtf_fn) as gtf:
        for i, record in enumerate(gtf):
            if record.startswith('#'):
                continue
            record = record.split('\t')
            chrom, _, feat_type, start, end, _, strand = record[:7]
            start = int(start) - 1
            end = int(end)
            if feat_type == 'exon':
                gene_id = get_gtf_attribute(record, 'gene_id')
                idx = (chrom, gene_id)
                if idx not in gtf_records:
                    gtf_records[idx] = []
                gtf_records[idx].append([start, end])
    gene_locs = {}
    for (chrom, gene_id), invs in  gtf_records.items():
        invs.sort()
        gene_locs[gene_id] = (chrom, (invs[0][0] + invs[-1][1]) // 2)
    return gene_locs


def flatten_intervals(invs):
    flattened = []
    inv_it = iter(invs)
    inv_start, inv_end = next(inv_it)
    for start, end in inv_it:
        if start <= inv_end:
            inv_end = max(inv_end, end)
        else:
            flattened.append([inv_start, inv_end])
            inv_start, inv_end = start, end
    if not flattened or flattened[-1] != [inv_start, inv_end]:
        flattened.append([inv_start, inv_end])
    return flattened


def parse_gtf(gtf_fn, quant_type='gene'):
    gtf_records = {}
    with open(gtf_fn) as gtf:
        for record in gtf:
            if record.startswith('#'):
                continue
            record = record.split('\t')
            chrom, _, feat_type, start, end, _, strand = record[:7]
            start = int(start) - 1
            end = int(end)
            if feat_type == 'exon':
                gene_id = get_gtf_attribute(record, 'gene_id')
                idx = (chrom, strand, gene_id)
                if idx not in gtf_records:
                    gtf_records[idx] = []
                gtf_records[idx].append([start, end])
    gene_ids = []
    gene_locs = []
    for (chrom, strand, gene_id), invs in gtf_records.items():
        invs.sort()
        if quant_type == 'gene':
            invs = flatten_intervals(invs)
        elif quant_type == 'gene_full':
            invs = [[invs[0][0], invs[-1][1]]]
        else:
            raise ValueError('unrecognised quant_type')
        gene_ids.append(gene_id)
        gene_locs.append((chrom, strand, invs))
    return gene_ids, gene_locs


def parse_query(query_str):
    try:
        chrom, *pos = re.search('^([\w\d]+):(\d+)(?:-(\d+))?$', query_str).groups()
        if pos[1]:
            pos = (int(pos[0]) + int(pos[1])) // 2
        else:
            pos = int(pos[0])
        assert pos > 0
    except:
        raise ValueError(f'Failed to parse haplotype query string "{query_str}"')
    return chrom, pos
