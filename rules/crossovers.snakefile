def get_co_input(wc):
    ref_name = config['reference_genotype']
    bams = expand(
        'aligned_data/diploid/{{cond}}.{parent1}__vs__{parent2}.sorted.bam',
        parent1=ref_name,
        parent2=[parent2 for parent2 in config['genome_fasta_fns'] if parent2 != ref_name]
    )
    return {'bams': bams, 'geno': 'cb_genotypes/{cond}.geno.tsv'}


rule detect_crossovers:
    input:
        unpack(get_co_input)
    output:
        model='cb_crossovers/{cond}.model.json',
        co='cb_crossovers/{cond}.tsv'
    threads: 5
    resources:
        mem_mb=10_000,
        queue='ioheavy'
    conda:
        'env_yamls/pomegranate.yaml'
    shell:
        '''
        python ../scripts/detect_crossovers.py \
          -g {input.geno} \
          --output-model-fn {output.model} \
          -o {output.co} \
          {input.bams}
        '''


rule get_cb_stats:
    input:
        geno='cb_genotypes/{cond}.geno.tsv',
        co='cb_crossovers/{cond}.tsv',
    output:
        stats='cb_stats/{cond}.tsv'
    threads: 5
    resources:
        mem_mb=10_000,
        queue='ioheavy'
    conda:
        'env_yamls/pomegranate.yaml'
    shell:
        '''
        python ../scripts/cb_stats.py \
          -g {input.geno} \
          -b {input.co} \
          -o {output.stats}
        '''


rule eqtl_analysis:
    input:
        solo_output='aligned_data/reference/{cond}_starsolo',
        co='cb_crossovers/{cond}.tsv',
        stats='cb_stats/{cond}.tsv'
    output:
        eqtls='eqtls/{cond}.tsv'
    threads: 30
    resources:
        mem_mb=lambda wc, threads: threads * 1024,
        queue='ioheavy'
    conda:
        'env_yamls/pomegranate.yaml'
    shell:
        '''
        python ../scripts/eqtl_analysis.py \
          -s {input.solo_output} \
          -b {input.co} \
          -w {input.stats} \
          -o {output.eqtls}
        '''


'''
def get_eqtl_input(wc):
    ref_name = config['reference_genotype']
    solo_output = expand(
        'aligned_data/diploid/{{cond}}.{parent1}__vs__{parent2}_starsolo',
        parent1=ref_name,
        parent2=[parent2 for parent2 in config['genome_fasta_fns'] if parent2 != ref_name]
    )
    return {'solo_output': solo_output,
            'co': 'cb_crossovers/{cond}.tsv',
            'stats': 'cb_stats/{cond}.tsv'}


def get_solo_params_eqtl(wc, input):
    genos = glob_wildcards('aligned_data/diploid/{cond}.{parent1}__vs__{parent2}_starsolo', input.solo_output).parent2
    param = ''
    for geno, fn in zip(genos, input.solo_output):
        param += f'-s {geno} {fn} '
    return param


rule eqtl_analysis:
    input:
        unpack(get_eqtl_input)
    output:
        eqtls='eqtls/{cond}.tsv'
    threads: 30
    resources:
        mem_mb=lambda wc, threads: threads * 1024,
        queue='ioheavy'
    conda:
        'env_yamls/pomegranate.yaml'
    params:
        solo_output=get_solo_params_eqtl
    shell:
        \'''
        python ../scripts/eqtl_analysis.py \
          {params.solo_output} \
          -b {input.co} \
          -w {input.stats} \
          -o {output.eqtls}
        \'''
'''