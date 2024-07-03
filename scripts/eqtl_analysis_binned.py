import os
import re
from copy import copy
import warnings
import itertools as it
from collections import defaultdict
from functools import partial

os.environ['OPENBLAS_NUM_THREADS'] = '1'

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score

from joblib import Parallel, delayed
import click

from singlecell_snptools.io import parse_gene_positions, parse_query
from singlecell_snptools.exprs import read_mmformat, lognormalise_exprs

LOD = 2 * np.log(10)
ORGANELLAR_CONTIGS = set(['ChrM', 'ChrC'])


def estimate_effective_haplotypes(hap_probs):
    '''
    Estimate the effective number of haplotypes in order to perform FWER correction
    of p values from eQTL tests (which are not independent due to LD)
    Uses the Li and Hi method (Heredity volume 95, pages221–227 (2005))
    '''
    H_eff = 0

    for chrom, h in hap_probs.T.groupby(level='chrom'):
        corr = h.T.corr(method='pearson')
        eig = np.abs(np.linalg.eigvalsh(corr))
        chrom_H_eff = np.sum((eig >= 1).astype(int) + (eig - np.floor(eig)))
        H_eff += round(chrom_H_eff)
    return H_eff


def find_pca_covariates(gene_exprs, haplotype_blocks, H_eff,
                        var_explained_threshold=0.01,
                        max_pcs=50,
                        fdr_threshold=0.05,
                        verbose=True):
    '''
    Finds principal components of gene_exprs that explain up to min_var_explained of total variance,
    removing any PCs that correlate with specific haplotypes
    (i.e. that are capturing differences in expression  caused by genetic differences)
    '''
    pca = PCA(n_components=max_pcs, svd_solver='arpack')
    principal_components = pca.fit_transform(gene_exprs.T)
    stop_idx = max_pcs - np.searchsorted(pca.explained_variance_ratio_[::-1], var_explained_threshold)
    if verbose:
        print(f'Using first {stop_idx} PCs that explain {pca.explained_variance_ratio_[:stop_idx].sum() * 100:.1f}% of gene expression variance')
    principal_components = pd.DataFrame(
        principal_components[:, :stop_idx],
        columns=[f'PC{i}' for i in range(1, stop_idx + 1)],
        index=gene_exprs.columns
    )

    pc_haplo_pvals = []
    pc_haplo_rsq = []
    for (chrom, pos), haplo in haplotype_blocks.items():
        pos_pc_pvals = []
        pos_pc_rsq = []
        for pc_idx, pc in principal_components.items():
            m, c, r, p, *_ = stats.linregress(pc, haplo)
            p = np.minimum(p * stop_idx * H_eff, 1.0) # bonferroni correction using H_eff
            pos_pc_pvals.append(p)
            pos_pc_rsq.append(r ** 2)
        pc_haplo_pvals.append(pos_pc_pvals)
        pc_haplo_rsq.append(pos_pc_rsq)
    pc_haplo_pvals = pd.DataFrame(pc_haplo_pvals, index=haplotype_blocks.columns, columns=principal_components.columns)
    pc_haplo_rsq = pd.DataFrame(pc_haplo_rsq, index=haplotype_blocks.columns, columns=principal_components.columns)
    pc_filt = np.logical_and(pc_haplo_pvals < fdr_threshold, pc_haplo_rsq >= var_explained_threshold).any(axis=0)
    if verbose:
        filt_pcs = ', '.join(principal_components.columns[pc_filt.values])
        print(f'Discarding PCs: [{filt_pcs},] which are correlated with haplotypes')
        var_explained = pca.explained_variance_ratio_[:stop_idx][~pc_filt.values].sum() * 100
        print(f'Remaining {stop_idx - pc_filt.sum()} PCs explain {var_explained:.2f}% of variance - these will be used as covariates')
    pc_haplo_stats = pc_haplo_pvals.join(pc_haplo_rsq, lsuffix='_pval', rsuffix='_r2')
    return principal_components.loc[:, ~pc_filt], principal_components.loc[:, pc_filt], pc_haplo_stats


def create_cell_type_clusters(principal_components, n_clusters='auto', covariance_type='spherical', random_state=None, verbose=True):
    max_s = 0
    if n_clusters == 'auto':
        for n in range(2, 10):
            gmm = GaussianMixture(
                n_components=n,
                covariance_type=covariance_type,
                random_state=random_state
            ).fit(principal_components)
            s = silhouette_score(principal_components, gmm.predict(principal_components))
            if s > max_s:
                max_s = s
                best_gmm = gmm
    else:
        best_gmm = GaussianMixture(
            n_components=n_clusters,
            random_state=random_state,
            covariance_type=covariance_type,
        ).fit(principal_components)
        max_s = silhouette_score(principal_components, best_gmm.predict(principal_components))
    if verbose:
        print(f'Identified {best_gmm.n_components} putative cell-type clusters with silhouette score {max_s:.2f}')
    celltype_clusters = pd.DataFrame(
        best_gmm.predict_proba(principal_components),
        index=principal_components.index,
        columns=[f'cluster{i}' for i in range(1, best_gmm.n_components + 1)]
    )
    return celltype_clusters


def fdr_fwer_correction(pos_eqtl_res, H_eff):
    pval_col_names = [col for col in pos_eqtl_res.columns if col.endswith('_pval')]
    for p_col in pval_col_names:
        # replace zero p values with smallest possible float64
        pos_eqtl_res[p_col] = np.maximum(pos_eqtl_res[p_col], np.finfo(np.float64).tiny)
        # correct gene-wise using BH, position-wise using effective number of haplotypes
        pos_eqtl_res[p_col] = np.minimum(
            fdrcorrection(pos_eqtl_res[p_col])[1] * H_eff, 1.0
        )
    return pos_eqtl_res


def likelihood_ratio(lln, llf):
    return -2 * (lln - llf)


def compare_lr_test(m_full, m_nest):
    
    Ln, Lf = m_nest.llf, m_full.llf
    dfn, dff = m_nest.df_model, m_full.df_model
    
    chi2_stat = likelihood_ratio(Ln, Lf)
    p = stats.chi2.sf(chi2_stat, dff - dfn)

    return chi2_stat, p


def haplotype_linear_model(gene_exprs, haplotype_blocks, cb_genotypes,
                           covariates, celltype_clusters,
                           control_haplotype, control_cis_haplotype, gene_locs, 
                           haplo_names, model_type='ols'):
    res = []
    n_cov = covariates.shape[1]
    haplo_pval_colnames = [f'{h}_pval' for h in haplo_names]
    haplo_coef_colnames = [f'{h}_coef' for h in haplo_names]
    lm = sm.Logit if model_type == 'logit' else partial(sm.OLS, hasconst=True)

    if cb_genotypes is not None:
        geno_dummies = pd.get_dummies(cb_genotypes.loc[:, ['geno']], drop_first=False).astype('float')
        geno_dummies.columns = geno_dummies.columns.str[5:]
        geno_dummies = geno_dummies[haplo_names]
        # drop first geno for covariates, keep for haplotypes
        all_covariates = pd.concat([geno_dummies.iloc[:, 1:], covariates], axis=1)
    else:
        geno_dummies = None
        all_covariates = covariates

    if control_haplotype is not None:
        if cb_genotypes is not None:
            control_haplo_dummies = geno_dummies.mul(control_haplotype, axis=0)
        else:
            control_haplo_dummies = control_haplotype.to_frame()
        all_covariates = pd.concat([control_haplo_dummies, covariates], axis=1)

    all_covariates.columns = all_covariates.columns + '_covar'

    if control_cis_haplotype:
        genes_to_test_for_haplotype, gene_cis_haplotypes = get_cis_haplotypes_and_linkage(
            gene_exprs.index.tolist(),
            gene_locs,
            haplotype_blocks,
            geno_dummies
        )
    else:
        genes_to_test_for_haplotype = {
            haplo: gene_exprs.index.tolist() for haplo in haplotype_blocks.columns
        }
        gene_cis_haplotype = None
    
    # generate nested models without haplotypes for likelihood ratio test
    nested_models = {}
    all_covariates_gene_specific = {}
    for gene_id, endog in gene_exprs.iterrows():
        if control_cis_haplotype:
            cis_haplo = gene_cis_haplotypes[gene_id]
            all_covariates_gene = pd.concat([all_covariates, cis_haplo], axis=1)
        else:
            all_covariates_gene = all_covariates

        m = lm(endog, sm.add_constant(all_covariates_gene)).fit(disp=0)

        all_covariates_gene_specific[gene_id] = all_covariates_gene
        nested_models[gene_id] = m

    if celltype_clusters is not None:
        haplo_x_celltype_names = [
            '_x_'.join(ixn) for ixn in it.product(haplo_names, celltype_clusters.columns)
        ]
        haplo_x_celltype_pval_colnames = [f'{h}_pval' for h in haplo_x_celltype_names]
        haplo_x_celltype_coef_colnames = [f'{h}_coef' for h in haplo_x_celltype_names]
        # generate nested models without haplotypes for likelihood ratio test
        nested_models_cluster = {}
        all_covariates_cluster_gene_specific = {}
        for gene_id, endog in gene_exprs.iterrows():
            if control_cis_haplotype:
                cis_haplo = gene_cis_haplotypes[gene_id]
                cis_haplo_x_cluster = pd.concat({f'{hap}_{gene_id}_cis_haplo':
                                                 celltype_clusters.mul(cis_haplo[f'{hap}_{gene_id}_cis_haplo'], axis=0)
                                                 for hap in haplo_names}, axis=1)
                cis_haplo_x_cluster.columns = ['_x_'.join(col).strip() for col in cis_haplo_x_cluster.columns.values]
                all_covariates_cluster = pd.concat([all_covariates, cis_haplo_x_cluster], axis=1)
            else:
                all_covariates_cluster = all_covariates
            all_covariates_cluster = pd.concat([celltype_clusters.iloc[:, 1:], all_covariates_cluster], axis=1)
            m = lm(endog, sm.add_constant(all_covariates_cluster)).fit(disp=0)
            all_covariates_cluster_gene_specific[gene_id] = all_covariates_cluster
            nested_models_cluster[gene_id] = m

    
    for (chrom, pos), haplo in haplotype_blocks.items():
        if geno_dummies is not None:
            haplo_dummies = geno_dummies.mul(haplo, axis=0)
        else:
            haplo_dummies = haplo.to_frame()
            haplo_dummies.columns = haplo_names
        
        for gene_id in genes_to_test_for_haplotype[(chrom, pos)]:
            endog = gene_exprs.loc[gene_id]
            exog = sm.add_constant(pd.concat([haplo_dummies, all_covariates_gene_specific[gene_id]], axis=1))
            gene_exprs = gene_exprs.loc[:, exog.index]

            if celltype_clusters is not None:
                haplo_x_celltype = pd.concat({hap: celltype_clusters.mul(haplo_dummies[hap], axis=0) for hap in haplo_names}, axis=1)
                haplo_x_celltype.columns = ['_x_'.join(col).strip() for col in haplo_x_celltype.columns.values]
                exog_cluster = sm.add_constant(pd.concat([haplo_x_celltype, all_covariates_cluster_gene_specific[gene_id]], axis=1))

            m_full = lm(endog, exog).fit(disp=0)
            m_nest = nested_models[gene_id]
            haplo_pvals = m_full.pvalues.reindex(haplo_names, fill_value=1.0).values
            haplo_coefs = m_full.params.reindex(haplo_names, fill_value=0.0).values
            lr, lrtest_pval = compare_lr_test(m_full, m_nest)
            
            pos_gene_res = [chrom, pos, gene_id, lr / LOD, lrtest_pval, *haplo_pvals, *haplo_coefs]
            
            if celltype_clusters is not None:
                m_cluster = lm(endog, exog_cluster).fit(disp=0)
                m_nest_cluster = nested_models_cluster[gene_id]
                lr_cluster, lrtest_cluster_pval = compare_lr_test(m_full, m_nest)
                haplo_x_cluster_pvals = m_cluster.pvalues.reindex(haplo_x_celltype_names, fill_value=1.0).values
                haplo_x_cluster_coefs = m_cluster.params.reindex(haplo_x_celltype_names, fill_value=0.0).values
                pos_gene_res += [lr_cluster / LOD, lrtest_cluster_pval, *haplo_x_cluster_pvals, *haplo_x_cluster_coefs]
            res.append(pos_gene_res)

    columns = [
        'chrom', 'pos', 'gene_id',
        'lod_score', 'lrt_pval',
        *haplo_pval_colnames, *haplo_coef_colnames
    ]
    if celltype_clusters is not None:
        columns += ['cluster_lod_score', 'cluster_lrt_pval', *haplo_x_celltype_pval_colnames, *haplo_x_celltype_coef_colnames]

    res = pd.DataFrame(res, columns=columns)
    return res


def parallel_haplo_lm(exprs_matrix, hap_probs, cb_genotypes, covariates, celltype_clusters,
                      control_haplotype, control_cis_haplotype, gene_locs, H_eff,
                      model_type='ols', processes=1):
    if cb_genotypes is not None:
        geno_names = sorted(cb_genotypes.geno.unique())
    else:
        geno_names = ['hap2',]

    args = (hap_probs, cb_genotypes, covariates, celltype_clusters,
            control_haplotype, control_cis_haplotype,
            gene_locs, geno_names, model_type)
    with Parallel(processes) as pool:
        res = pool(
            delayed(haplotype_linear_model)(exprs_chunk, *args)
            for exprs_chunk in np.array_split(exprs_matrix, processes)
        )
    res = pd.concat(res)
    res = res.groupby(['chrom', 'pos']).apply(fdr_fwer_correction, H_eff=H_eff)
    return res


def find_closest_haplotype(hap_probs, chrom, pos):
    hp = hap_probs[chrom]
    return hp.iloc[:, np.argmin(np.abs(hp.columns.astype(int) - pos))]


def find_closely_linked_haplotypes(hap_probs, test_haplo, min_r2=0.95):
    r2 = hap_probs.apply(lambda hap: stats.pearsonr(hap, test_haplo)[0] ** 2, axis=0)
    return r2 >= min_r2


def get_haplotype_covariates(hap_probs, control_haplotype, remove_r2=0.95):
    chrom, pos = parse_query(control_haplotype)
    print(f'Adding haplotype at {chrom}:{pos:d} as a model covariate')
    control_haplotype = find_closest_haplotype(hap_probs, chrom, pos).rename('control_haplo')
    r2_mask = find_closely_linked_haplotypes(hap_probs, control_haplotype, remove_r2)
    print(f'Removing {sum(r2_mask)} haplotype bins that are in linkage with {chrom}:{pos:d}')
    return control_haplotype, hap_probs.loc[:, ~r2_mask]


def get_cis_haplotypes_and_linkage(gene_ids, gene_locs, hap_probs, geno_dummies=None, min_r2=0.95):
    genes_to_test_for_haplotype = defaultdict(list)
    gene_cis_haplotype = {}
    for gene_id in gene_ids:
        chrom, pos = gene_locs[gene_id]
        cis_haplo = find_closest_haplotype(hap_probs, chrom, pos)
        linked_haplotypes = find_closely_linked_haplotypes(hap_probs, cis_haplo, min_r2)
        for (chrom, pos) in linked_haplotypes.index[~linked_haplotypes]:
            genes_to_test_for_haplotype[(chrom, pos)].append(gene_id)
        if geno_dummies is not None:
            cis_haplo = geno_dummies.mul(cis_haplo, axis=0)
            cis_haplo.columns = geno_dummies.columns + f'{gene_id}_cis_haplo'
        else:
            cis_haplo = cis_haplo.rename(f'hap2_{gene_id}_cis_haplo').to_frame()
        gene_cis_haplotype[gene_id] = cis_haplo
    return genes_to_test_for_haplotype, gene_cis_haplotype


def parse_peak_id(peak_id):
    chrom, range_ = peak_id.split(':')
    start, end = range_.split('-')
    return chrom, (int(start) + int(end)) // 2


@click.command()
@click.option('-s', '--starsolo-dir', required=True, multiple=True)
@click.option('-p', '--hap-probs-fn', required=True, multiple=True)
@click.option('-w', '--whitelist-fn', required=True)
@click.option('-o', '--output-prefix', required=True)
@click.option('--min-cells-exprs', default=0.05)
@click.option('--celltype-x-haplotype-interaction/--no-interaction', default=False)
@click.option('--clusters', required=False, default='auto')
@click.option('--control-haplotype', required=False, default=None)
@click.option('--control-cis-haplotype/--no-control-cis-haplotype', default=False)
@click.option('--gtf-fn', required=False, default=None)
@click.option('--control-haplotype-r2', required=False, default=0.95)
@click.option('--data-type', default='rna', type=click.Choice(['rna', 'atac'], case_sensitive=False))
@click.option('-g', '--genotype-filter', default=None)
@click.option('-n', '--processes', required=False, default=1)
@click.option('-r', '--random-state', required=False, default=101)
def main(starsolo_dir, hap_probs_fn, whitelist_fn, output_prefix, min_cells_exprs,
         celltype_x_haplotype_interaction, clusters,
         control_haplotype, control_cis_haplotype, gtf_fn, control_haplotype_r2,
         data_type, genotype_filter, processes, random_state):
    print('Reading cb stats')
    try:
        cb_stats = pd.read_csv(
            whitelist_fn,
            header=0,
            index_col=0,
            sep='\t',
            usecols=['cb', 'geno', 'cluster'],
        )
    except ValueError:
        cb_stats = pd.read_csv(
            whitelist_fn,
            header=0,
            index_col=0,
            sep='\t',
            usecols=['cb', 'cluster'],
        )
    if genotype_filter is not None:
        cb_stats = cb_stats.query('geno == @genotype_filter')
        assert len(cb_stats), f'No cell barcodes left after filtering for genotype {genotype_filter}'
        cb_stats.drop(columns=['geno'], inplace=True)
    cb_whitelist = sorted(cb_stats.query('cluster == "hq"').index.tolist())
    print(f'{len(cb_whitelist)} high quality cell barcodes identified & whitelisted')

    print('Reading expression matrix')
    exprs_mat = read_mmformat(starsolo_dir, cb_whitelist)
    exprs_mat = exprs_mat.loc[:, cb_whitelist]
    if data_type == 'rna':
        ln_exprs_mat, _ = lognormalise_exprs(exprs_mat)
    elif data_type == 'atac':
        # binarise 
        ln_exprs_mat = (exprs_mat > 0).astype(int)
    else:
        raise NotImplementedError('Only RNA and ATAC supported')
    
    hap_probs = pd.concat((
        pd.read_csv(fn, sep='\t', header=[0, 1], index_col=[0, 1])
        for fn in hap_probs_fn
    ), axis=0)
    hap_probs = hap_probs[
        [col for col in hap_probs.columns.get_level_values(0).unique() if col not in ORGANELLAR_CONTIGS]
    ]
    hap_probs.columns.set_names(['chrom', 'pos'], inplace=True)
    hap_probs.index = hap_probs.index.droplevel(0)
    hap_probs = hap_probs.loc[cb_whitelist]

    if control_haplotype is not None:
        control_haplotype, hap_probs = get_haplotype_covariates(
            hap_probs, control_haplotype, remove_r2=control_haplotype_r2
        )

    if control_cis_haplotype:
        if gtf_fn is None and data_type == 'rna':
            raise ValueError('GTF file must be provided to use --control-cis-haplotype')
        if data_type == 'rna':
            gene_positions = parse_gene_positions(gtf_fn)
        else:
            gene_positions = {peak_id: parse_peak_id(peak_id) for peak_id in ln_exprs_mat.index.tolist()}
    else:
        gene_positions = None

    H_eff = estimate_effective_haplotypes(hap_probs)
    print(f'Effective haplotype number is {H_eff} - this value will be used for FWER correction')

    try:
        cb_genotypes = cb_stats.loc[cb_whitelist, ['geno']]
    except KeyError:
        cb_genotypes = None

    (principal_components, 
     principal_components_discarded,
     pc_haplo_stats) = find_pca_covariates(ln_exprs_mat, hap_probs, H_eff)
    if celltype_x_haplotype_interaction:
        if len(principal_components.columns):
            celltype_clusters = create_cell_type_clusters(
                principal_components,
                clusters,
                random_state=random_state
            )
        else:
            warnings.warn('No principal components that do not correlate with haplotypes. Cell types cannot be clustered. Possibly you only have one cell type.')
            celltype_clusters = None
    else:
        celltype_clusters = None

    if data_type == 'rna':
        # filter out genes which are not expressed in enough cells to be worth testing
        ln_exprs_mat = ln_exprs_mat[(exprs_mat > 0).mean(axis=1) > min_cells_exprs]
    elif data_type == 'atac':
        # filter out peaks where either too few or nearly all cells have fragments
        ln_exprs_mat = ln_exprs_mat[exprs_mat.mean(axis=1) > min_cells_exprs]
        ln_exprs_mat = ln_exprs_mat[exprs_mat.mean(axis=1) < (1 - min_cells_exprs)]
    print(f'{len(ln_exprs_mat)} {"genes" if data_type == "rna" else "peaks"} will be tested')
    
    res = parallel_haplo_lm(
        ln_exprs_mat, hap_probs,
        cb_genotypes,
        principal_components,
        celltype_clusters,
        control_haplotype,
        control_cis_haplotype,
        gene_positions,
        H_eff=H_eff,
        model_type='ols' if data_type == 'rna' else 'logit',
        processes=processes
    )
    res.to_csv(
        f'{output_prefix}.eqtls.tsv',
        sep='\t',
        index=False,
        header=True,
        float_format='%.4g'
    )
    covars = principal_components.join(celltype_clusters) if celltype_clusters is not None else principal_components
    covars.to_csv(
        f'{output_prefix}.covars.tsv',
        sep='\t',
        index=True,
        header=True,
        float_format='%.4g'
    )
    principal_components_discarded.to_csv(
        f'{output_prefix}.discarded_pcs.tsv',
        sep='\t',
        index=True,
        header=True,
        float_format='%.4g'
    )
    pc_haplo_stats.to_csv(
        f'{output_prefix}.pc_haplo_stats.tsv',
        sep='\t',
        index=True,
        header=True,
        float_format='%.4g'
    )


if __name__ == '__main__':
    main()