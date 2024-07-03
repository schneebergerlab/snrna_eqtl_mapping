def get_nb_input(wc):
    ref_name = config['datasets'][wc.cond]['reference_genotype']
    parent2_accs = config['datasets'][wc.cond]['parent2_accessions']
    bams = expand(
        'aligned_data/diploid/{{cond}}.{parent1}__vs__{parent2}.sorted.bam',
        parent1=ref_name,
        parent2=parent2_accs
    )
    cos = expand(
        'cb_crossovers/{{cond}}.{parent1}__vs__{parent2}.co_probs.tsv',
        parent1=ref_name,
        parent2=parent2_accs
    )
    exprs = expand(
        'aligned_data/diploid/{{cond}}.{parent1}__vs__{parent2}_genecounts',
        parent1=ref_name,
        parent2=parent2_accs
    )
    input_ = {
        'bams': bams,
        'ref_fasta': ancient(config['genome_fasta_fns'][ref_name]),
        'parent2_fastas': [ancient(config['genome_fasta_fns'][acc]) for acc in parent2_accs],
        'gtf': ancient(config['gtf_fns'][ref_name]),
        'cb_stats': 'cb_stats/{cond}.tsv',
        'crossovers': cos,
        'eqtls': 'eqtls/{cond}.all.eqtls.tsv',
        'exprs': exprs,
    }
    if len(parent2_accs) > 1:
        input_['vcf'] = 'annotation/diploid/vcf/{cond}.cellsnplite.vcf.gz'
    return input_


rule expression_figures:
    '''
    Perform expression analysis and generate more permissive CB whitelist
    '''
    input:
        unpack(get_nb_input)
    output:
        figures=directory('figures/expression_figures_{cond}'),
    conda:
        'env_yamls/notebooks.yaml'
    notebook:
        'notebook_templates/expression_figures.{wildcards.cond}.py.ipynb'


rule genotyping_figures:
    '''
    Perform expression analysis and generate more permissive CB whitelist
    '''
    input:
        unpack(get_nb_input)        
    output:
        figures=directory('figures/genotyping_figures_{cond}'),
    conda:
        'env_yamls/notebooks.yaml'
    notebook:
        'notebook_templates/genotyping_figures.{wildcards.cond}.py.ipynb'


rule crossover_figures:
    '''
    Perform expression analysis and generate more permissive CB whitelist
    '''
    input:
        unpack(get_nb_input)
    output:
        figures=directory('figures/crossover_figures_{cond}'),
    threads: 5
    conda:
        'env_yamls/notebooks.yaml'
    notebook:
        'notebook_templates/crossover_figures.{wildcards.cond}.py.ipynb'


rule eqtl_figures:
    '''
    Perform expression analysis and generate more permissive CB whitelist
    '''
    input:
        unpack(get_nb_input),
    output:
        figures=directory('figures/eqtl_figures_{cond}'),
    conda:
        'env_yamls/notebooks.yaml'
    notebook:
        'notebook_templates/eqtl_figures.{wildcards.cond}.py.ipynb'


rule hotspot_figures:
    '''
    Perform expression analysis and generate more permissive CB whitelist
    '''
    input:
        unpack(get_nb_input)
    output:
        figures=directory('figures/hotspot_figures_{cond}'),
    conda:
        'env_yamls/notebooks.yaml'
    notebook:
        'notebook_templates/hotspot_figures.{wildcards.cond}.py.ipynb'
      