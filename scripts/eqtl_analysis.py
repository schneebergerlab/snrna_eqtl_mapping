import os
import re
from copy import copy

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
from sklearn.decomposition import PCA

from joblib import Parallel, delayed
import click

from singlecell_snptools.exprs import read_mmformat
from singlecell_snptools.crossovers import get_haplotype_blocks


LOD = 2 * np.log(10)


def haplotype_linear_model(gene_exprs, haplotype_blocks, covariates, haplo_names, ref_name='col0'):
    res = []
    n_cov = covariates.shape[1]
    wald_expr = ' , '.join(haplo_names)
    haplo_pval_colnames = [f'{h}_pval' for h in haplo_names]
    # generate nested models without haplotypes for likelihood ratio test
    nested_models = {}
    for gene_id, endog in gene_exprs.iterrows():
        m = sm.OLS(endog, sm.add_constant(covariates), hasconst=True).fit(disp=0)
        nested_models[gene_id] = m
    for (chrom, start, end), haplo in haplotype_blocks.iterrows():
        haplo_dummies = (
            pd.get_dummies(haplo, drop_first=False)
                .astype(float)
                .drop(ref_name, axis=1)
        )
        exog = sm.add_constant(pd.concat([haplo_dummies, covariates], axis=1))
        gene_exprs = gene_exprs.loc[:, exog.index]
        for gene_id, endog in gene_exprs.iterrows():
            m_full = sm.OLS(endog, exog, hasconst=True).fit(disp=0)
            m_nest = nested_models[gene_id]
            wald_pval = float(m_full.wald_test(wald_expr, scalar=True).pvalue)
            haplo_pvals = m_full.pvalues.reindex(haplo_names, fill_value=1.0).values
            lr, lrtest_pval, _ = m_full.compare_lr_test(m_nest)
            res.append([chrom, start, end, gene_id, lr / LOD, lrtest_pval, wald_pval, *haplo_pvals])
    res = pd.DataFrame(res, columns=['chrom', 'start', 'end', 'gene_id',
                                     'lod_score', 'lrt_pval', 'wald_pval',
                                     *haplo_pval_colnames])
    return res


def parallel_haplo_lm(exprs_matrix, co_invs, cb_genotypes, covariates, ref_name='col0', processes=1):
    haplotype_blocks = get_haplotype_blocks(co_invs)
    geno_names = sorted(cb_genotypes.geno.unique())
    geno_dummies = pd.get_dummies(cb_genotypes.loc[:, ['geno']], drop_first=True).astype('float')
    all_covariates = pd.concat([geno_dummies, covariates], axis=1)
    with Parallel(processes) as pool:
        res = pool(
            delayed(haplotype_linear_model)(exprs_chunk, haplotype_blocks, all_covariates, geno_names, ref_name)
            for exprs_chunk in np.array_split(exprs_matrix, processes)
        )
    res = pd.concat(res)
    pval_col_names = [col for col in res.columns if col.endswith('_pval')]
    for p_col in pval_col_names:
        f_col = p_col.rsplit('_', 1)[0] + '_fdr'
        _, res[f_col] = fdrcorrection(res[p_col])
    res = res.sort_values(['chrom', 'start', 'gene_id'])
    return res


@click.command()
@click.option('-s', '--starsolo-dir', required=True)
@click.option('-b', '--co-bed-fn', required=True)
@click.option('-w', '--whitelist-fn', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('--pca-covariates', required=False, default=2)
@click.option('--min-cells-exprs', default=0.1)
@click.option('--ref-accession', required=False, default='col0')
@click.option('-p', '--processes', required=False, default=1)
def main(starsolo_dir, co_bed_fn, whitelist_fn, output_fn, pca_covariates, min_cells_exprs, ref_accession, processes):
    exprs_mat = read_mmformat(
        os.path.join(starsolo_dir, 'Gene/raw/matrix.mtx'),
        os.path.join(starsolo_dir, 'Gene/raw/barcodes.tsv'),
        os.path.join(starsolo_dir, 'Gene/raw/features.tsv'),
    )
    co_invs = pd.read_csv(
        co_bed_fn,
        sep='[\t|]',
        engine='python',
        names=['chrom', 'start', 'end', 'cb', 'haplo'],
        usecols=[0, 1, 2, 3, 4],
        dtype={'chrom': 'str'}
    )
    cb_stats = pd.read_csv(
        whitelist_fn,
        header=0,
        index_col=0,
        sep='\t',
        usecols=['cb', 'geno', 'cluster'],
    )
    cb_whitelist = sorted(cb_stats.query('cluster == "hq"').index.tolist())
    cb_genotypes = cb_stats.loc[cb_whitelist, ['geno']]
    
    # filter out cell barcodes that failed genotyping
    exprs_mat = exprs_mat.loc[:, cb_whitelist]
    co_invs = co_invs[co_invs.cb.isin(cb_whitelist)]
    # filter out genes which are not expressed in any cells
    exprs_mat = exprs_mat[(exprs_mat > 0).mean(axis=1) > min_cells_exprs]

    pca = PCA(n_components=pca_covariates)
    principal_components = pd.DataFrame(
        pca.fit_transform(exprs_mat.T),
        columns=[f'PC{i}' for i in range(1, pca_covariates + 1)],
        index=exprs_mat.columns
    )

    res = parallel_haplo_lm(
        exprs_mat, co_invs, 
        cb_genotypes,
        principal_components,
        ref_name=ref_accession,
        processes=processes
    )
    res.to_csv(output_fn, sep='\t', index=False, header=True, float_format='%.4g')


if __name__ == '__main__':
    main()