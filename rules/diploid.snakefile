def get_star_diploid_index_input(wc):
    parent1, parent2 = wc.comp.split('__vs__')
    return {
        'fasta_fn': ancient(config['genome_fasta_fns'][parent1]),
        'gtf_fn': ancient(config['gtf_fns'][parent1]),
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
        STAR \
          --runThreadN {threads} \
          --runMode genomeGenerate \
          --genomeDir {output} \
          --genomeFastaFiles {input.fasta_fn} \
          --sjdbGTFfile {input.gtf_fn} \
          --genomeTransformVCF {input.vcf_fn} \
          --genomeTransformType Diploid \
          --genomeSAindexNbases 12 \
          --sjdbOverhang {params.overhang}
        '''


rule STAR_diploid:
    '''
    map reads with STAR spliced aligner
    '''
    input:
        read_barcode='demuxed_data/{cond}.{comp}.1_barcode.fastq.gz',
        mate='demuxed_data/{cond}.{comp}.2.fastq.gz',
        index='annotation/diploid/{comp}_STAR_index',
        barcode_whitelist=ancient(config['barcode_whitelist']),
    output:
        bam='aligned_data/diploid/{cond}.{comp}.sorted.bam',
        bai='aligned_data/diploid/{cond}.{comp}.sorted.bam.bai',
        stats='aligned_data/diploid/{cond}.{comp}.sorted.bamstats',
        sjdb='aligned_data/diploid/{cond}.{comp}.sjdb.tsv',
        solo_output=directory('aligned_data/diploid/{cond}.{comp}_starsolo')
    params:
        sort_mem=lambda wc, resources: int((resources.mem_mb - 4096) * 1e6),
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
          --readFilesIn "$TOPDIR/{input.mate}" "$TOPDIR/{input.read_barcode}" \
          --readFilesCommand "zcat" \
          --soloType "CB_UMI_Simple" \
          --soloCBwhitelist "$TOPDIR/{input.barcode_whitelist}" \
          --soloUMIlen 12 \
          --soloUMIdedup "1MM_Directional_UMItools" \
          --outFilterMultimapNmax 2 \
          --outFilterIntronMotifs RemoveNoncanonical \
          --alignSJoverhangMin 12 \
          --alignSJDBoverhangMin 4 \
          --outFilterMismatchNmax 2 \
          --alignIntronMin 60 \
          --alignIntronMax 20000 \
          --outSAMtype BAM SortedByCoordinate \
          --outBAMsortingBinsN 150 \
          --limitBAMsortRAM {params.sort_mem} \
          --outSAMattributes NH HI AS nM CB UB ha \
          --genomeTransformOutput SAM SJ Quant

        cd $TOPDIR
        mv ${{STAR_TMP_DIR}}/Aligned.sortedByCoord.out.bam {output.bam}
        mv ${{STAR_TMP_DIR}}/Log.progress.out {log.progress}
        mv ${{STAR_TMP_DIR}}/Log.final.out {log.final}
        mv ${{STAR_TMP_DIR}}/Log.out {log.main}
        mv ${{STAR_TMP_DIR}}/SJ.out.tab {output.sjdb}
        mv ${{STAR_TMP_DIR}}/Solo.out {output.solo_output}
          
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.stats}
        rm -rf $STAR_TMP_DIR
        '''
