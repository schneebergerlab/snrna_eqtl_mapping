import re
from collections import defaultdict
from statistics import median_low

import numpy as np
import pandas as pd
from scipy.io import mmread

import pysam
from ncls import NCLS


def get_chrom_sizes(vcf_fn):
    with pysam.VariantFile(vcf_fn) as vcf:
        chrom_sizes = {k: v.length for k, v in vcf.header.contigs.items()}
    return chrom_sizes


def get_chrom_sizes_bam(bam_fn):
    with pysam.AlignmentFile(bam_fn) as bam:
        chrom_sizes = {k: bam.get_reference_length(k) for k in bam.references}
    return chrom_sizes


def parse_snps(vcf_fn, minq=10):
    with pysam.VariantFile(vcf_fn) as vcf:
        accessions = list(vcf.header.samples)
        chrom_sizes = {k: v.length for k, v in vcf.header.contigs.items()}
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
        return parsed_snps, accessions, chrom_sizes


def read_cb_whitelist(barcode_fn):
    with open(barcode_fn) as f:
        cb_whitelist = [cb.strip() for cb in f.readlines()]
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
    snps, accessions, chrom_sizes = parse_snps(vcf_fn)
    snps = pd.MultiIndex.from_tuples(snps, names=['chrom', 'pos', 'ref', 'alt', 'ref_accs', 'alt_accs'])
    alt_mm = pd.DataFrame.sparse.from_spmatrix(alt_mm, columns=barcodes, index=snps)
    ref_mm = pd.DataFrame.sparse.from_spmatrix(ref_mm, columns=barcodes, index=snps)
    counts = pd.concat({'ref_count': ref_mm, 'alt_count': alt_mm}, axis=1)
    counts.columns = counts.columns.reorder_levels([1, 0])
    return counts, accessions, chrom_sizes


def read_whitelist_itrees(syn_bed_fn):
    syn_regions = pd.read_csv(
        syn_bed_fn,
        sep='\t',
        names=['chrom', 'start', 'end', 'geno'],
        dtype={'chrom': str}
    )

    syn_region_itrees = {}
    for (geno, chrom), invs in syn_regions.groupby(['geno', 'chrom']):
        syn_region_itrees[(geno, chrom)] = NCLS(
            invs.start.values,
            invs.end.values,
            np.arange(len(invs))
        )
    return syn_region_itrees


def filter_snps_by_itree(snps, geno, whitelist_itrees):
    filtered_snps = []
    for chrom, chrom_snps in snps.groupby('chrom'):
        filt_idx, _ = whitelist_itrees[(geno, chrom)].all_overlaps_both(
            chrom_snps.pos.values - 1,
            chrom_snps.pos.values,
            np.arange(len(chrom_snps))
        )
        filt_idx = np.unique(filt_idx)
        filtered_snps.append(chrom_snps.iloc[filt_idx])
    return pd.concat(filtered_snps)


def collapse_short_marker_distances(chrom_co_markers, cluster_dist=250):
    chrom_co_markers = chrom_co_markers.sort_values(by='pos')
    cluster_ids = (chrom_co_markers.pos.diff().bfill() > cluster_dist).cumsum()
    return chrom_co_markers.groupby(cluster_ids, as_index=False).agg(
        {'chrom': 'first', 'pos': median_low, 'ref_count': 'sum', 'alt_count': 'sum'}
    )


def get_lr(co_markers, pseudo=0.1):
    return np.log(co_markers.alt_count + pseudo) - np.log(co_markers.ref_count + pseudo)


def convert_to_co_markers(cb_snp_counts, geno, whitelist_region_itrees=None, cluster_dist=250):
    m = cb_snp_counts.apply(lambda snp: geno in snp.alt_accs, axis=1)
    cb_co_markers = (cb_snp_counts.loc[m, ['chrom', 'pos', 'ref_count', 'alt_count']]
                                  .sort_values(['chrom', 'pos']))
    if whitelist_region_itrees is not None:
        cb_co_markers = filter_snps_by_itree(cb_co_markers, geno, whitelist_region_itrees)
    cb_co_markers = (cb_co_markers.groupby('chrom', as_index=False)
                                  .apply(collapse_short_marker_distances, cluster_dist=cluster_dist))
    cb_co_markers = cb_co_markers.reset_index(drop=True)
    cb_co_markers['lr'] = get_lr(cb_co_markers)
    return cb_co_markers


def iter_cb(snp_counts, genotypes=None, as_type='snps', cluster_dist=250, whitelist_region_itrees=None):
    if as_type == 'co' and genotypes is None:
        raise ValueError()
    cb_list = list(set(snp_counts.columns.get_level_values(0)))
    np.random.shuffle(cb_list)
    for cb in cb_list:
        cb_snp_counts = snp_counts[cb]
        cb_snp_counts = cb_snp_counts.iloc[
            np.unique(np.concatenate([
                cb_snp_counts.ref_count.values.sp_index.indices,
                cb_snp_counts.alt_count.values.sp_index.indices,
            ]))
        ]
        cb_snp_counts = cb_snp_counts.sparse.to_dense().reset_index()
        if as_type == 'snps':
            yield cb, None, cb_snp_counts
        elif as_type == 'co':
            cb_geno = genotypes[cb]
            cb_co_markers = convert_to_co_markers(
                cb_snp_counts,
                cb_geno,
                whitelist_region_itrees=whitelist_region_itrees,
                cluster_dist=cluster_dist
            )
            yield cb, cb_geno, cb_co_markers


def get_gtf_attribute(gtf_record, attribute):
    try:
        attr = re.search(f'{attribute} "(.+?)";', gtf_record[8]).group(1)
    except AttributeError:
        raise ValueError(
            f'Could not parse attribute {attribute} '
            f'from GTF with feature type {record[2]}'
        )
    return attr


def get_gene_locations(gtf_fn):
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
        gene_locs[gene_id] = (chrom, invs[0][0], invs[-1][1])
    return pd.DataFrame.from_dict(gene_locs, orient='index', columns=['chrom', 'start', 'end'])


def get_co_markers_from_star_diploid_bam(bam_fn, cluster_dist=250):
    co_markers = []
    seen_umi_cb = defaultdict(set)
    with pysam.AlignmentFile(bam_fn) as bam:
        for aln in bam.fetch():
            if aln.is_secondary or aln.is_supplementary or aln.is_duplicate:
                continue
            hap = aln.get_tag('ha')
            cb, umi = aln.get_tag('CB'), aln.get_tag('UB')
            if hap != 0:
                if umi not in seen_umi_cb[cb]:
                    co_markers.append([cb, aln.reference_name, aln.reference_start, hap])
                    seen_umi_cb[cb].add(umi)
        co_markers = pd.DataFrame(co_markers, columns=['cb', 'chrom', 'pos', 'hap'])
        co_markers = (co_markers.groupby(['cb', 'chrom', 'pos'])
                                .hap.value_counts()
                                .unstack(level=-1)
                                .fillna(0)
                                .astype(int)
                                .reset_index())
        co_markers.columns = ['cb', 'chrom', 'pos', 'ref_count', 'alt_count']
        co_markers = (co_markers.groupby(['cb', 'chrom'])
                                .apply(collapse_short_marker_distances, cluster_dist=cluster_dist)
                                .reset_index(level=['cb'])
                                .reset_index(drop=True))
        co_markers['lr'] = get_lr(co_markers)
    return co_markers