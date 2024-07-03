import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
import seaborn as sns


def chrom_axes_boilerplate(chrom_sizes, figsize=(18, 5), nrows=1, row_ratios=None, span_features=None, span_kwargs=None):

    if nrows == 'nchroms':
        fig, axes = plt.subplots(
            figsize=figsize,
            ncols=len(chrom_sizes),
            nrows=len(chrom_sizes),
            height_ratios=list(chrom_sizes.values()),
            width_ratios=list(chrom_sizes.values()),
            sharey='row',
            sharex='col'
        )
        for chrom, ax in zip(chrom_sizes, axes[:, 0]):
            ax.set_ylim(0, chrom_sizes[chrom])
            ax.set_yticks(np.arange(0, chrom_sizes[chrom], 1e7))
            ax.set_yticklabels([int(i // 1e6) for i in np.arange(0, chrom_sizes[chrom], 1e7)])
            ax.set_ylabel(f'Chromosome {chrom} (Mb)')
    else:
        fig, axes = plt.subplots(
            figsize=figsize,
            ncols=len(chrom_sizes),
            nrows=nrows,
            height_ratios=row_ratios,
            width_ratios=list(chrom_sizes.values()),
            sharey='row',
            sharex='col'
        )
    axes = axes.reshape(-1, len(chrom_sizes))
    for chrom, ax in zip(chrom_sizes, axes[-1]):
        ax.set_xlim(0, chrom_sizes[chrom])
        ax.set_xticks(np.arange(0, chrom_sizes[chrom], 1e7))
        ax.set_xticklabels([int(i // 1e6) for i in np.arange(0, chrom_sizes[chrom], 1e7)])
        ax.set_xlabel(f'Chromosome {chrom} (Mb)')

    if span_kwargs is None:
        span_kwargs = {'alpha': 0.3, 'color': '#252525'}

    for row in axes:
        for chrom, ax in zip(chrom_sizes, row):
            if span_features is not None:
                for s, e in span_features[chrom]:
                    ax.axvspan(s, e, **span_kwargs)

    plt.tight_layout()
    return fig, axes


def ref_alt_rugplot(co_markers, chrom_size, bin_size, ax=None, max_height=10):
    pos = np.arange(len(co_markers)) * bin_size
    rwgt = co_markers[:, 0].astype(np.int16)
    rwgt[rwgt > max_height] = max_height
    awgt = co_markers[:, 1].astype(np.int16)
    awgt[awgt > max_height] = max_height
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    ax.vlines(pos, np.repeat(0, len(pos)), rwgt, color='#0072b2', zorder=0)
    ax.vlines(pos, np.repeat(0, len(pos)), np.negative(awgt), color='#d55e00', zorder=0)
    ax.plot([0, chrom_size], [0, 0], ls='-', color='#252525')
    return ax


DEFAULT_CMAP = LinearSegmentedColormap.from_list('hap_cmap', ['#0072b2', '#d55e00'])
YLIM_OFFSET = 1
XLIM_OFFSET = 1e4


def single_cell_co_plot(co_markers, hap_probs, chrom_sizes, bin_size=25_000,
                        show_mesh_prob=True, show_line_prob=True, annotate_co_number=True,
                        max_height=10, cmap=DEFAULT_CMAP, figsize=(18, 4)):
    fig, axes = chrom_axes_boilerplate(chrom_sizes, figsize=figsize)
    axes = axes.squeeze()
    axes[0].set_ylabel('Marker coverage')
    for chrom, ax in zip(chrom_sizes, axes):
        ax.set_xlim(-XLIM_OFFSET, chrom_sizes[chrom] + XLIM_OFFSET)
        ax.set_ylim(-max_height - YLIM_OFFSET, max_height + YLIM_OFFSET)
        ref_alt_rugplot(co_markers[chrom], chrom_sizes[chrom], bin_size, ax=ax, max_height=max_height)
        hp = hap_probs.loc[chrom]
        x_pos = hp.index.values.astype(int)
        if show_mesh_prob:
            c = ax.pcolormesh(
                np.insert(x_pos, len(x_pos), chrom_sizes[chrom]),
                np.array([-max_height - YLIM_OFFSET, max_height + YLIM_OFFSET]),
                hp.values.reshape(1, -1),
                cmap=cmap,
                norm=Normalize(0, 1),
                alpha=0.5,
                zorder=-2,
                rasterized=True
            )
            if chrom == list(chrom_sizes)[-1]:
                cax = plt.colorbar(c, ax=axes)
                cax.ax.set_ylabel('Predicted haplotype')
                cax.ax.set_yticks([0, 0.5, 1])
                cax.ax.set_yticklabels(['Hap1', '50/50', 'Hap2'])
                cax.ax.set_rasterized(True)
        if show_line_prob:
            ax.plot(
                x_pos,
                np.negative(hp.values - 0.5) * (max_height * 2 - 1),
                color='#EEEEEE' if show_mesh_prob else '#252525',
                lw=2,
                zorder=-1
            )
        if annotate_co_number:
            n_co = np.abs(np.diff(hp.values)).sum()
            ax.annotate(text=f'{n_co.sum():.2f} COs', xy=(0.05, 0.1), xycoords='axes fraction')
    return fig, axes


def create_eqtl_plotter(eqtl_res, sig_eqtl, gene_locs, col_mapper, palette, lod_score=False):
    
    def _eqtl_plotter(gene_id, title=None):
        gene_eqtl_res = (eqtl_res.query('gene_id == @gene_id')
                                 .melt(id_vars=['chrom', 'pos'],
                                       value_vars=list(col_mapper),
                                       var_name='hue',
                                       value_name='y'))
        if not len(gene_eqtl_res):
            raise ValueError('Gene was not tested')
        gene_eqtl_res = gene_eqtl_res.assign(
            hue=gene_eqtl_res.hue.map(col_mapper),
            y=np.negative(np.log10(gene_eqtl_res.y)) if not lod_score else gene_eqtl_res.y
        )
        if gene_locs is not None:
            gene_chrom, gene_pos = gene_locs[gene_id]
        else:
            gene_chrom, gene_pos = None, None
        chrom_sizes = eqtl_res.groupby('chrom').pos.max().to_dict()
        fig, axes = chrom_axes_boilerplate(chrom_sizes, figsize=(12, 3.5))
        axes = axes.ravel()
        for (chrom, chrom_eqtl_res), ax in zip(gene_eqtl_res.groupby('chrom'), axes):
            sns.lineplot(
                x='pos',
                y='y',
                hue='hue',
                data=chrom_eqtl_res,
                drawstyle='steps-post',
                ax=ax,
                hue_order=list(col_mapper.values()),
                palette=palette
            )
            if chrom == gene_chrom:
                ax.axvline(gene_pos, ls='-', color='#252525', zorder=-1)
            for _, se in sig_eqtl.query('gene_id == @gene_id & chrom_eqtl == @chrom').iterrows():
                ax.axvline(se.pos_eqtl, ls='--', color='#777777', zorder=-1)
                ax.axvspan(se.start_eqtl, se.end_eqtl, color='#777777', zorder=-1, alpha=0.2)

        if not lod_score:
            axes[0].set_ylabel('Negative log10 FDR')
        else:
            axes[0].set_ylabel('LOD score')
        for ax in axes[:-1]:
            ax.legend_.remove()
        axes[-1].legend_.set_title('')

        if title is None:
            title = gene_id
        fig.suptitle(title)

        plt.tight_layout(pad=0.2)
        return axes
    return _eqtl_plotter


def _get_haplotype_position(eqtls, hap_probs, haplotype_type, chrom=None, pos=None, gene_chrom=None, gene_pos=None):
    if haplotype_type == 'cis':
        assert gene_chrom is not None and gene_pos is not None
        haplotype = eqtls.loc[eqtls.query('chrom == @gene_chrom').eval('abs(pos - @gene_pos)').idxmin()]
        name = 'Cis-haplotype'
    elif haplotype_type == 'best_hit':
        haplotype = eqtls.loc[eqtls.lod_score.idxmax()]
        name= f'Haplotype at {haplotype.chrom}:{haplotype.pos}'
    elif haplotype_type == 'specified':
        assert chrom is not None
        if pos is not None:
            haplotype = eqtls.loc[eqtls.query('chrom == @chrom').eval('abs(pos - @pos)').idxmin()]
        else:
            # take best hit on chrom
            haplotype = eqtls.loc[eqtls.query('chrom == @chrom').lrt_pval.idxmin()]
        name= f'Haplotype at {haplotype.chrom}:{haplotype.pos}'
    hap = hap_probs.loc[:, pd.IndexSlice[haplotype.chrom, str(haplotype.pos)]]
    return (hap > 0.5).map({False: 'Parent 1', True: 'Parent 2'}).rename(name)


def create_haplotype_expression_plotter(exprs_mat, eqtls, hap_probs, cb_stats, gene_locs):
    cb_whitelist = cb_stats.query('cluster == "hq"').index.values
    exprs_mat = exprs_mat.loc[:, cb_whitelist]
    hap_probs = hap_probs.loc[cb_whitelist]
    cb_stats = cb_stats.loc[cb_whitelist]

    def _hap_exprs_plot(gene_id, x='haplotype', hue=None, *,
                        haplotype='cis', haplotype_hue=None,
                        haplotype_chrom=None, haplotype_pos=None,
                        haplotype_hue_chrom=None, haplotype_hue_pos=None,
                        ax=None, figsize=(8, 5), palette=None, cb_subset_query=None,
                        **violin_kws):
        gene_eqtls = eqtls.query('gene_id == @gene_id')
        gene_chrom, gene_pos = gene_locs[gene_id]
        hap = _get_haplotype_position(
            gene_eqtls, hap_probs, haplotype,
            chrom=haplotype_chrom, pos=haplotype_pos,
            gene_chrom=gene_chrom, gene_pos=gene_pos
        )
        if haplotype_hue is not None:
            hue_hap = _get_haplotype_position(
                gene_eqtls, hap_probs, haplotype_hue,
                chrom=haplotype_hue_chrom, pos=haplotype_hue_pos,
                gene_chrom=gene_chrom, gene_pos=gene_pos
            )
        else:
            hue_hap = None

        y = exprs_mat.loc[gene_id]
        if x == 'haplotype':
            x = hap
        elif x == 'genotype':
            x = cb_stats['geno'].rename('Genotype')
        else:
            x = cb_stats[x].rename(x.capitalize())

        if hue is None:
            hue = None
        elif hue == 'haplotype':
            hue = hue_hap
        elif hue == 'genotype':
            hue = cb_stats['geno'].rename('Genotype')
        else:
            hue = cb_stats[hue].rename(hue.capitalize())

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if cb_subset_query is not None:
            mask = cb_stats.eval(cb_subset_query)
            y = y[mask]
            x = x[mask]
            hue = hue[mask]
            
        sns.violinplot(
            y=y,
            x=x,
            hue=hue,
            ax=ax,
            **violin_kws
        )
        ax.set_ylabel(f'{gene_id} log2 gene expression')

        return ax
    return _hap_exprs_plot