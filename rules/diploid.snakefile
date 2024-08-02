def get_star_diploid_index_input(wc):
    parent1, parent2 = wc.comp.split('__vs__')
    return {
        'fasta_fn': ancient(config['genome_fasta_fns'][parent1]),
        'organellar_fasta_fn': ancient(config['organellar_genome_fasta_fn']),
        'gtf_fn': ancient(config['gtf_fns'][parent1]),
        'organellar_gtf_fn': ancient(config['organellar_gtf_fn']),
        'vcf_fn': 'annotation/diploid/vcf/{comp}.snvs.vcf'
    }


rule build_STAR_index_diploid:
    '''Create the index required for alignment with STAR'''
    input:
        unpack(get_star_diploid_index_input)
    output:
        directory('annotation/diploid/{comp}_STAR_index')
    threads: 12
    resources:
        mem_mb=12 * 1024,
        queue='ioheavy'
    params:
        overhang=150,
    conda:
        'env_yamls/star.yaml'
    shell:
        '''
        mkdir {output};
        cat {input.fasta_fn} {input.organellar_fasta_fn} > \
          annotation/diploid/{wildcards.comp}_fullref.fa
        cat {input.gtf_fn} {input.organellar_gtf_fn} > \
          annotation/diploid/{wildcards.comp}_fullref.gtf
        STAR \
          --runThreadN {threads} \
          --runMode genomeGenerate \
          --genomeDir {output} \
          --genomeFastaFiles annotation/diploid/{wildcards.comp}_fullref.fa \
          --sjdbGTFfile annotation/diploid/{wildcards.comp}_fullref.gtf \
          --genomeTransformVCF {input.vcf_fn} \
          --genomeTransformType Diploid \
          --genomeSAindexNbases 12 \
          --sjdbOverhang {params.overhang}

        rm annotation/diploid/{wildcards.comp}_fullref.fa \
           annotation/diploid/{wildcards.comp}_fullref.gtf
        '''


def sample_name_subset(cond):
    sample_names = glob_wildcards(
        'raw_data/{sample_name}.1.fastq.gz'
    ).sample_name
    cond_sample_names = [sn for sn in sample_names if sn.rsplit('_', 1)[0] == cond]
    return sorted(cond_sample_names)


def expand_sample_name_from_cond(pattern):
    def _expand_sn(wc):
        return expand(
            pattern,
            sample_name=sample_name_subset(wc.cond)
        )
    return _expand_sn


def STAR_diploid_input(wc):
    parent2_accs = config['datasets'][wc.cond]['parent2_accessions']
    if len(parent2_accs) > 1:
        return dict(
            read_barcode=['demuxed_data/{cond}.{comp}.1_barcode.fastq.gz'],
            mate=['demuxed_data/{cond}.{comp}.2.fastq.gz'],
            index='annotation/diploid/{comp}_STAR_index',
            barcode_whitelist=ancient(config['barcode_whitelist']),
        )
    else:
        return dict(
            read_barcode=expand_sample_name_from_cond('trimmed_data/{sample_name}.1_barcode.fastq.gz')(wc),
            mate=expand_sample_name_from_cond('raw_data/{sample_name}.2.fastq.gz')(wc),
            index='annotation/diploid/{comp}_STAR_index',
            barcode_whitelist=ancient(config['barcode_whitelist']),
        )


rule STAR_diploid:
    '''
    map reads with STAR spliced aligner
    '''
    input:
        unpack(STAR_diploid_input)
    output:
        bam='aligned_data/diploid/{cond}.{comp}.sorted.bam',
        bai='aligned_data/diploid/{cond}.{comp}.sorted.bam.bai',
        stats='aligned_data/diploid/{cond}.{comp}.sorted.bamstats',
        sjdb='aligned_data/diploid/{cond}.{comp}.sjdb.tsv',
    params:
        read_barcode=lambda wc, input: ','.join(f'${{TOPDIR}}/{fn}' for fn in input.read_barcode),
        mate=lambda wc, input: ','.join(f'${{TOPDIR}}/{fn}' for fn in input.mate),
        sort_mem=lambda wc, resources: (resources.mem_mb - 4096) * 1_000_000,
        n_files=lambda wc, threads: threads * 150 + 200
    log:
        progress='logs/{cond}.{comp}.STAR_progress.log',
        final='logs/{cond}.{comp}.STAR_final.log',
        main='logs/{cond}.{comp}.STAR.log'
    threads: 24
    resources:
        mem_mb=lambda wildcards, threads: (threads + 4) * 2048,
        queue='ioheavy'
    conda:
        'env_yamls/star.yaml'
    shell:
        '''
        TOPDIR=$(pwd)
        STAR_TMP_DIR="aligned_data/diploid/{wildcards.cond}.{wildcards.comp}.tmpdir"
        mkdir -p $STAR_TMP_DIR
        cd $STAR_TMP_DIR
        ulimit -n {params.n_files}
        STAR \
          --runThreadN {threads} \
          --genomeDir "$TOPDIR/{input.index}" \
          --readFilesIn "{params.mate}" "{params.read_barcode}" \
          --readFilesCommand "zcat" \
          --soloType "CB_UMI_Simple" \
          --soloCBwhitelist "$TOPDIR/{input.barcode_whitelist}" \
          --soloUMIlen 12 \
          --soloUMIdedup "NoDedup" \
          --outFilterMultimapNmax 2 \
          --outFilterIntronMotifs RemoveNoncanonical \
          --alignSJoverhangMin 12 \
          --alignSJDBoverhangMin 4 \
          --outFilterMismatchNmax 4 \
          --alignIntronMin 60 \
          --alignIntronMax 20000 \
          --outSAMtype BAM SortedByCoordinate \
          --outBAMsortingBinsN 150 \
          --limitBAMsortRAM {params.sort_mem} \
          --outSAMattributes NH HI AS nM NM CB UB UR ha \
          --genomeTransformOutput SAM SJ Quant

        cd $TOPDIR
        mv ${{STAR_TMP_DIR}}/Aligned.sortedByCoord.out.bam {output.bam}
        mv ${{STAR_TMP_DIR}}/Log.progress.out {log.progress}
        mv ${{STAR_TMP_DIR}}/Log.final.out {log.final}
        mv ${{STAR_TMP_DIR}}/Log.out {log.main}
        mv ${{STAR_TMP_DIR}}/SJ.out.tab {output.sjdb}
          
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.stats}
        rm -rf $STAR_TMP_DIR
        '''


rule filter_star_reads_not_overlapping_variants:
    input:
        bam='aligned_data/diploid/{cond}.{comp}.sorted.bam',
        vcf='annotation/diploid/vcf/{comp}.snvs.vcf',
    output:
        bam='aligned_data/diploid/{cond}.{comp}.filt.bam',
        bai='aligned_data/diploid/{cond}.{comp}.filt.bam.bai',
    resources:
        mem_mb=10_000,
        queue='ioheavy'
    conda:
        'env_yamls/samtools.yaml'
    shell:
        '''
        awk -v OFS='\t' \
          '$0 !~ /^#/ {{print $1, $2 - 2, $2 + length($4), $3}}' \
          {input.vcf} > {output.bam}.snvs.bed

        samtools view -b \
          -q 255 \
          -e '[ha] != 0 && [nM] == 0 && [CB] !~ "-"' \
          -L {output.bam}.snvs.bed {input.bam} \
          > {output.bam}

        samtools index {output.bam}
        rm {output.bam}.snvs.bed
        '''


def get_gene_count_input(wc):
    parent1, parent2 = wc.comp.split('__vs__')
    return {
        'bam_fn': 'aligned_data/diploid/{cond}.{comp}.sorted.bam',
        'gtf_fn': ancient(config['gtf_fns'][parent1]),
        'organellar_gtf_fn': ancient(config['organellar_gtf_fn']),
        'whitelist_fn': 'aligned_data/reference/{cond}.cb_whitelist.txt'
    }


rule get_diploid_gene_counts:
    input:
        unpack(get_gene_count_input)
    output:
        directory('aligned_data/diploid/{cond}.{comp}_genecounts')
    resources:
        mem_mb=25_000,
        queue='ioheavy'
    conda:
        'env_yamls/pysam.yaml'
    threads: 25
    shell:
        '''
        python ../scripts/sc_gene_counts.py -n {threads} \
          -g <(cat {input.gtf_fn} {input.organellar_gtf_fn}) \
          -w {input.whitelist_fn} \
          -b {input.bam_fn} \
          -o {output}
        '''