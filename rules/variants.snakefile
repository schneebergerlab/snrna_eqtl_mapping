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
         minimap2 --eqx -t 8 -ax asm20 \
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
    resources:
        mem_mb=lambda wc, threads: 2048 * threads,
        queue='ioheavy'
    conda:
        'env_yamls/syri.yaml'
    shell:
        '''
        syri -F B -f --hdrseq --dir annotation/diploid/vcf --prefix {wildcards.comp}. \
          -c {input.bam} -q {input.parent2} -r {input.parent1}
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


def get_merge_input(wc):
    ref_name = config['reference_genotype']
    vcfs = expand(
        'annotation/diploid/vcf/{parent1}__vs__{parent2}.syri.vcf',
        parent1=ref_name,
        parent2=[parent2 for parent2 in config['genome_fasta_fns'] if parent2 != ref_name]
    )
    return {'vcfs': vcfs, 'fai': ancient(config['genome_fasta_fns'][config['reference_genotype']] + '.fai')}


rule merge_snps_for_genotyping:
    input:
        unpack(get_merge_input)
    output:
        'annotation/reference/vcf/{cond}.syri.vcf.gz',
    params:
        gz_vcf=lambda wc, input: [fn + '.gz' for fn in input.vcfs],
        gz_vcf_idx=lambda wc, input: [fn + '.gz.tbi' for fn in input.vcfs],
        filt_vcf=lambda wc, input: [fn + '.filt.vcf.gz' for fn in input.vcfs],
        filt_vcf_idx=lambda wc, input: [fn + '.filt.vcf.gz.tbi' for fn in input.vcfs],
    conda:
        'env_yamls/pysam.yaml'
    resources:
        mem_mb=10_000
    shell:
        '''
        awk -v OFS='\\t' '{{print $1, "0", $2}}' {input.fai} > {output}.syn.bed
        for VCF in {input.vcfs}; do
          awk -v OFS='\\t' \
            'substr($3, 1, 5) == "SYNAL" \
            {{split($8, INFO, /[;=]/); print $1, $2, INFO[2]}}' $VCF | \
          bedtools merge -i stdin | \
          bedtools intersect -a {output}.syn.bed -b stdin > {output}.syn.tmp.bed
          mv {output}.syn.tmp.bed {output}.syn.bed
        done
        for VCF in {input.vcfs}; do
          BASENAME="${{VCF##*__vs__}}"
          PARENT2="${{BASENAME%%.syri.vcf}}"
          bgzip -c $VCF > "${{VCF}}.gz"
          tabix -p vcf "${{VCF}}.gz"
          bcftools view -R {output}.syn.bed -i 'TYPE="snp"' "${{VCF}}.gz" | \
          awk -v SAMPLE=$PARENT2 -v OFS='\\t' \
            '{{if ($0 ~ /^##/) {{print}} \
               else if ($0 ~ /^#/) {{print "##FORMAT=<ID=GT,Number=1,Type=String>"; print $0, "FORMAT", SAMPLE}} \
               else {{print $1, $2, $3, $4, $5, $6, $7, ".", "GT", "1/1"}}}}' | \
          bgzip > "${{VCF}}.filt.vcf.gz"
          tabix -p vcf "${{VCF}}.filt.vcf.gz"
        done
        rm {output}.syn.bed
        bcftools merge -0 -O z {params.filt_vcf} > {output}.temp.vcf.gz
        rm {params.gz_vcf}
        rm {params.gz_vcf_idx}
        rm {params.filt_vcf}
        rm {params.filt_vcf_idx}
        tabix -p vcf {output}.temp.vcf.gz
        bcftools view -M2 -O z {output}.temp.vcf.gz > {output}
        rm {output}.temp.vcf.gz*
        tabix -p vcf {output}
        '''