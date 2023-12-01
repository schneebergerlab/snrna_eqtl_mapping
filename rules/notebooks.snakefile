rule expression_figures:
    '''
    Perform expression analysis and generate more permissive CB whitelist
    '''
    input:
        solo_output='aligned_data/reference/{cond}_starsolo',
        cb_stats='cb_stats/{cond}.tsv'
    output:
        figures=directory('figures/expression_figures_{cond}'),
    conda:
        'env_yamls/notebooks.yaml'
    notebook:
        'notebook_templates/expression_figures.py.ipynb'


rule genotyping_figures:
    '''
    Perform expression analysis and generate more permissive CB whitelist
    '''
    input:
        cb_stats='cb_stats/{cond}.tsv',
        vcf='annotation/reference/vcf/{cond}.syri.vcf.gz'
    output:
        figures=directory('figures/genotyping_figures_{cond}'),
    conda:
        'env_yamls/notebooks.yaml'
    notebook:
        'notebook_templates/genotyping_figures.py.ipynb'


def get_co_nb_input(wc):
    ref_name = config['reference_genotype']
    bams = expand(
        'aligned_data/diploid/{{cond}}.{parent1}__vs__{parent2}.sorted.bam',
        parent1=ref_name,
        parent2=[parent2 for parent2 in config['genome_fasta_fns'] if parent2 != ref_name]
    )
    return {
        'bams': bams,
        'crossovers': 'cb_crossovers/{cond}.tsv',
        'cb_stats': 'cb_stats/{cond}.tsv',
        'vcf': 'annotation/reference/vcf/{cond}.syri.vcf.gz',
        'gtf_fn': ancient(config['gtf_fns'][config['reference_genotype']])
    }


rule crossover_figures:
    '''
    Perform expression analysis and generate more permissive CB whitelist
    '''
    input:
        unpack(get_co_nb_input)
    output:
        figures=directory('figures/crossover_figures_{cond}'),
    threads: 5
    conda:
        'env_yamls/notebooks.yaml'
    notebook:
        'notebook_templates/crossover_figures.py.ipynb'


rule eqtl_figures:
    '''
    Perform expression analysis and generate more permissive CB whitelist
    '''
    input:
        eqtls='eqtls/{cond}.tsv',
        cb_stats='cb_stats/{cond}.tsv',
        vcf='annotation/reference/vcf/{cond}.syri.vcf.gz',
        gtf_fn=ancient(config['gtf_fns'][config['reference_genotype']]),
    output:
        figures=directory('figures/eqtl_figures_{cond}'),
    conda:
        'env_yamls/notebooks.yaml'
    notebook:
        'notebook_templates/eqtl_figures.py.ipynb'
      