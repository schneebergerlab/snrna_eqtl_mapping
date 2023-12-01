import numpy as np
import pandas as pd
import pomegranate as pm

import click


def read_cb_stats(geno_fn, co_fn):
    geno = pd.read_csv(
        geno_fn,
        sep='\t',
        header=0,
        names=['cb', 'geno', 'nmarkers', 'geno_logprob'],
        index_col='cb'
    )
    co_invs = pd.read_csv(
        co_fn,
        sep='[\t|]', engine='python',
        dtype={'chrom': str},
        names=['chrom', 'start', 'end', 'cb', 'geno', 'nmarkers', 'marker_range', 'hmm_logprob', 'strand']
    )
    co_invs['min_hmm_logprob'] = co_invs.hmm_logprob / co_invs.nmarkers
    co_invs['marker_density'] = co_invs.nmarkers / (co_invs.end - co_invs.start)
    co_stats = co_invs.groupby('cb').agg(dict(
        nmarkers='sum',
        hmm_logprob='sum',
        min_hmm_logprob='min',
        marker_density=lambda m: m.std() / m.mean(),
        chrom=lambda m: (m.value_counts() - 1).sum()
    ))
    co_stats['mean_hmm_logprob'] = co_stats.hmm_logprob / co_stats.nmarkers
    geno['n_crossovers'] = geno.index.map(co_stats.chrom)
    geno['co_score_mean'] = geno.index.map(co_stats.mean_hmm_logprob)
    geno['co_score_min'] = geno.index.map(co_stats.min_hmm_logprob)
    geno['density_cov'] = geno.index.map(co_stats.marker_density)
    geno['nmarkers'] = np.log10(geno.nmarkers)
    return geno


USESTATS = [
    'geno_logprob', 'co_score_mean', 'nmarkers',
]

def cluster_hq_barcodes(barcode_stats, init_geno_thresh=1.5, stat_names=USESTATS):
    barcode_stats = barcode_stats.copy()[stat_names]
    # good cells should in general have high geno probs
    init_label = (barcode_stats.geno_logprob > init_geno_thresh)
    dists = []
    for i in (0, 1):
        X_i = barcode_stats.loc[init_label == i]
        dists.append(pm.IndependentComponentsDistribution([
            pm.NormalDistribution(X_i[sn].mean(), X_i[sn].std())
            for sn in stat_names
        ]))
    gmm = pm.GeneralMixtureModel(dists)
    gmm.fit(barcode_stats.values, inertia=0.9)
    pred = gmm.predict(barcode_stats.values)
    return pd.Series(pred, index=barcode_stats.index).map({0: 'artefact', 1: 'hq'})


@click.command()
@click.option('-g', '--geno-fn', required=True)
@click.option('-b', '--co-bed-fn', required=True)
@click.option('-o', '--output-fn', required=True)
def main(geno_fn, co_bed_fn, output_fn):
    cb_stats = read_cb_stats(geno_fn, co_bed_fn)
    cb_stats['cluster'] = cluster_hq_barcodes(cb_stats)
    cb_stats.to_csv(output_fn, sep='\t')


if __name__ == '__main__':
    main()