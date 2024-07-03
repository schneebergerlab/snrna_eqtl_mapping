def get_comps(wc):
    parent2_accs = config['datasets'][wc.cond]['parent2_accessions']
    ref_name = config['datasets'][wc.cond]['reference_genotype']
    return expand('{p1}__vs__{p2}', p1=ref_name, p2=parent2_accs)


rule train_rhmm:
    input:
        bam=lambda wc: expand(
            'aligned_data/diploid/{cond}.{comp}.filt.bam',
            cond=wc.cond,
            comp=get_comps(wc)
        ),
        cb_whitelist='aligned_data/reference/{cond}.cb_whitelist.txt'
    output:
        model='cb_crossovers/{cond}.all_samples.model.json',
    threads: 20
    resources:
        mem_mb=50_000,
        queue='ioheavy'
    conda:
        'env_yamls/pomegranate.yaml'
    shell:
        '''
        python ../scripts/detect_crossovers_binned.py --train-only -p {threads} \
          -w {input.cb_whitelist} \
          {input.bam} \
          --output-model-fn {output.model}
        '''


rule detect_crossovers:
    input:
        bam='aligned_data/diploid/{cond}.{comp}.filt.bam',
        model='cb_crossovers/{cond}.all_samples.model.json',
        cb_whitelist='aligned_data/reference/{cond}.cb_whitelist.txt'
    output:
        cb_stats='cb_crossovers/{cond}.{comp}.cb_stats.tsv',
        co_probs='cb_crossovers/{cond}.{comp}.co_probs.tsv',
    threads: 20
    resources:
        mem_mb=50_000,
        queue='ioheavy'
    conda:
        'env_yamls/pomegranate.yaml'
    shell:
        '''
        python ../scripts/detect_crossovers_binned.py -p {threads} \
          -w {input.cb_whitelist} \
          -m {input.model} \
          -o cb_crossovers/{wildcards.cond}.{wildcards.comp} \
          {input.bam}
        '''

def cb_stats_input(wc):
    parent2_accs = config['datasets'][wc.cond]['parent2_accessions']
    if len(parent2_accs) > 1:
        return dict(
            geno='cb_genotypes/{cond}.geno.tsv',
            cb_stats=expand('cb_crossovers/{{cond}}.{comp}.cb_stats.tsv', comp=get_comps(wc)),
            co_probs=expand('cb_crossovers/{{cond}}.{comp}.co_probs.tsv', comp=get_comps(wc)),
        )
    else:
        return dict(
            cb_stats=expand('cb_crossovers/{{cond}}.{comp}.cb_stats.tsv', comp=get_comps(wc)),
            co_probs=expand('cb_crossovers/{{cond}}.{comp}.co_probs.tsv', comp=get_comps(wc)),
        )


rule get_cb_stats:
    input:
        unpack(cb_stats_input)
    output:
        stats='cb_stats/{cond}.tsv'
    threads: 5
    resources:
        mem_mb=10_000,
        queue='ioheavy'
    conda:
        'env_yamls/notebooks.yaml'
    params:
        stats_flag=lambda wc, input: ' -s '.join(input.cb_stats),
        probs_flag=lambda wc, input: ' -p '.join(input.co_probs),
        geno_flag=lambda wc, input: f'-g {input.geno}' if hasattr(input, 'geno') else ''
    notebook:
        'notebook_templates/cb_stats.{wildcards.cond}.py.ipynb'


def get_eqtl_analysis_input(wc):
    if wc.comp == 'all':
        comps = get_comps(wc)
        return {
            'solo_output': expand('aligned_data/diploid/{{cond}}.{comp}_genecounts', comp=comps),
            'co_probs': expand('cb_crossovers/{{cond}}.{comp}.co_probs.tsv', comp=comps),
            'stats': 'cb_stats/{cond}.tsv'
        }
    else:
        return {
            'solo_output': ['aligned_data/diploid/{cond}.{comp}_genecounts'],
            'co_probs': ['cb_crossovers/{cond}.{comp}.co_probs.tsv'],
            'stats': 'cb_stats/{cond}.tsv'
        }


rule eqtl_analysis:
    input:
        unpack(get_eqtl_analysis_input)
    output:
        eqtls='eqtls/{cond}.{comp}.eqtls.tsv'
    threads: 30
    resources:
        mem_mb=lambda wc, threads: threads * 2048,
        queue='ioheavy'
    conda:
        'env_yamls/pomegranate.yaml'
    params:
        solo_flag=lambda wc, input: ' -s '.join(input.solo_output),
        probs_flag=lambda wc, input: ' -p '.join(input.co_probs),
        geno_flag=lambda wc: f'-g {wc.comp.split("__vs__")[1]}' if wc.comp != 'all' else ''
    shell:
        '''
        python ../scripts/eqtl_analysis_binned.py -n {threads} \
          --celltype-x-haplotype-interaction \
          {params.geno_flag} \
          -s {params.solo_flag} \
          -p {params.probs_flag} \
          -w {input.stats} \
          -o eqtls/{wildcards.cond}.{wildcards.comp}
        '''
