def get_wga_input(wc):
    parent1, parent2 = wc.comp.split('__vs__')
    return {
        'parent1': ancient(config['genome_fasta_fns'][parent1]),
        'parent2': ancient(config['genome_fasta_fns'][parent2]),
    }


rule minimap2_wga:
    input:
        unpack(get_wga_input)
    output:
        'annotation/diploid/wga/{comp}.bam'
    threads:
        12
    resources:
        mem_mb=lambda wc, threads: 1024 * threads,
        queue='ioheavy'
    conda:
        'env_yamls/minimap2.yaml'
    shell:
        '''
         minimap2 --eqx -t 8 -ax asm20 -z1000,100 --no-end-flt \
           {input.parent1} {input.parent2}  | \
        samtools view -bS | \
        samtools sort -o - - > {output}
        samtools index {output}
        '''


def get_syri_input(wc):
    input_ = get_wga_input(wc)
    input_['bam'] = 'annotation/diploid/wga/{comp}.bam'
    return input_


rule run_syri:
    input:
        unpack(get_syri_input)
    output:
        vcf='annotation/diploid/vcf/{comp}.syri.vcf',
        syri='annotation/diploid/vcf/{comp}.syri.out',
    resources:
        mem_mb=lambda wc, threads: 2048 * threads,
        queue='ioheavy'
    conda:
        'env_yamls/msyd.yaml'
    params:
        acc=lambda wc: wc.comp.split('__vs__')[1]
    shell:
        '''
        syri -F B -f --hdrseq --dir annotation/diploid/vcf --prefix {wildcards.comp}. \
          -c {input.bam} -q {input.parent2} -r {input.parent1} --samplename {params.acc}
        '''


def get_msyd_input(wc):
    accs = config['datasets'][wc.cond]['parent2_accessions']
    ref_name = config['datasets'][wc.cond]['reference_genotype']
    return {
        'bams': expand('annotation/diploid/wga/{ref_name}__vs__{accession}.bam',
                       ref_name=ref_name, accession=accs),
        'syri': expand('annotation/diploid/vcf/{ref_name}__vs__{accession}.syri.out',
                       ref_name=ref_name, accession=accs),
        'vcf': expand('annotation/diploid/vcf/{ref_name}__vs__{accession}.syri.vcf',
                       ref_name=ref_name, accession=accs),
        'fasta': [config['genome_fasta_fns'][acc] for acc in accs],
    }
    

rule msyd_input:
    input:    
        unpack(get_msyd_input)
    output:
        cfg=temp('annotation/diploid/vcf/{cond}.msyd_config.tsv')
    group: 'msyd'
    run:
        ref_name = config['datasets'][wildcards.cond]['reference_genotype']
        with open(output.cfg, 'w') as f:
            f.write('#name\taln\tsyri\tvcf\tgenome\n')
            for accession in config['datasets'][wildcards.cond]['parent2_accessions']:
                f.write(
                    f'{accession}\t'
                    f'annotation/diploid/wga/{ref_name}__vs__{accession}.bam\t'
                    f'annotation/diploid/vcf/{ref_name}__vs__{accession}.syri.out\t'
                    f'annotation/diploid/vcf/{ref_name}__vs__{accession}.syri.vcf\t'
                    f'{config["genome_fasta_fns"][accession]}\n'
                )


rule run_msyd:
    input:
        cfg='annotation/diploid/vcf/{cond}.msyd_config.tsv',
        ref=lambda wc: ancient(config['genome_fasta_fns'][config['datasets'][wc.cond]['reference_genotype']])
    output:
        pff='annotation/diploid/vcf/{cond,\w+}.pansyn.pff',
        vcf='annotation/diploid/vcf/{cond,\w+}.vcf',
    group: 'msyd'
    conda:
        'env_yamls/msyd.yaml'
    shell:
        '''
        msyd call --core -i {input.cfg} -r {input.ref} -o {output.pff} -m {output.vcf}
        '''


rule convert_vcf_for_star:
    input:
        'annotation/diploid/vcf/{comp}.syri.vcf'
    output:
        'annotation/diploid/vcf/{comp}.snvs.vcf'
    conda:
        'env_yamls/pysam.yaml'
    resources:
        mem_mb=30_000
    shell:
        '''
        python ../scripts/convert_syri_output.py {input} {output}
        '''


rule filter_snps_for_star_consensus:
    input:
        vcf='annotation/diploid/vcf/{cond}.vcf'
    output:
        vcf='annotation/diploid/vcf/{cond,\w+}.consensus.vcf',
    conda:
        'env_yamls/pysam.yaml'
    params:
        min_ac=lambda wc: len(config['datasets'][wc.cond]['parent2_accessions']) // 2 + 1,
        max_indel_size=50,
    shell:
        r'''
        bcftools annotate \
          --exclude 'ALT ~ "CORESYN"' \
          --remove "FORMAT/CHR,FORMAT/START,FORMAT/END,INFO/PID" \
          {input.vcf} |
        grep -v "^##ALT" | grep -v "^##INFO" | \
        bcftools sort | \
        bcftools filter -S0 -e 'GT=="."' | \
        bcftools view -G -M2 \
          --min-ac {params.min_ac} \
          -e "STRLEN(REF)>{params.max_indel_size} || STRLEN(ALT)>{params.max_indel_size}" \
        > {output.vcf}
        '''


rule filter_snps_for_cellsnplite:
    input:
        vcf='annotation/diploid/vcf/{cond}.vcf'
    output:
        vcf='annotation/diploid/vcf/{cond,\w+}.cellsnplite.vcf.gz',
    conda:
        'env_yamls/pysam.yaml'
    shell:
        r'''
        bcftools annotate \
          --exclude 'ALT ~ "CORESYN"' \
          --remove "FORMAT/CHR,FORMAT/START,FORMAT/END,INFO/PID" \
          {input.vcf} |
        grep -v "^##ALT" | grep -v "^##INFO" | \
        bcftools sort | \
        bcftools filter -S0 -e 'GT=="."' | \
        bcftools view -M2 \
          -i 'TYPE="snp"' \
          -O z \
        > {output.vcf}
        tabix -p vcf {output.vcf}
        '''