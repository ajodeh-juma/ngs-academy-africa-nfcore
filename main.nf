#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// define workflow parameters
params.reads = "$projectDir/data/ggal/*_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.outdir = "results"

log.info """\
        R N A S E Q - N F   P I P E L I N E    
        ===================================
        transcriptome: ${params.transcriptome_file}
        reads        : ${params.reads}
        outdir       : ${params.outdir}
        """
        .stripIndent()




workflow {

    // data channel for transcriptome file
    ch_transcriptome = file(params.transcriptome_file, checkIfExists: true)
    
    // data channel for the input reads
    ch_raw_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // run fastqc process
    FASTQC(ch_raw_reads)
    ch_fastqc_multiqc = FASTQC.out.zip

    // run fastp process
    FASTP(ch_raw_reads)
    ch_preprocessed_reads = FASTP.out.reads
    ch_fastp_multiqc = FASTP.out.json
    
    // run salmon processes
    INDEX(ch_transcriptome)
    QUANT(INDEX.out.index, ch_preprocessed_reads)
    ch_salmon_multiqc = QUANT.out.results

    // run multiqc process
    MULTIQC( 
        ch_fastqc_multiqc.collect{it[1]}.ifEmpty([]),
        ch_salmon_multiqc.collect{it[1]}.ifEmpty([]),
        ch_fastp_multiqc.collect{it[1]}.ifEmpty([])
    )
}

process FASTQC {

    tag "$sample_id"

    publishDir "${params.outdir}/fastqc", mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip

    script:
    def prefix   = "${sample_id}"
    // Add soft-links to original FastQs for consistent naming in pipeline
    """
    [ ! -f  ${prefix}_1.fq ] && ln -s ${reads[0]} ${prefix}_1.fq
    [ ! -f  ${prefix}_2.fq ] && ln -s ${reads[1]} ${prefix}_2.fq

    fastqc --threads $task.cpus ${prefix}_1.fq ${prefix}_2.fq
    """
}

process FASTP {

    tag "$sample_id"

    publishDir "${params.outdir}/fastp", mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*.trim.fq')      , emit: reads
    tuple val(sample_id), path('*.json')         , emit: json
    tuple val(sample_id), path('*.html')         , emit: html
    tuple val(sample_id), path('*.log')          , emit: log

    script:
    def prefix   = "${sample_id}"
    // Add soft-links to original FastQs for consistent naming in pipeline
    """
    [ ! -f  ${prefix}_1.fq ] && ln -s ${reads[0]} ${prefix}_1.fq
    [ ! -f  ${prefix}_2.fq ] && ln -s ${reads[1]} ${prefix}_2.fq

    fastp \\
        --in1 ${prefix}_1.fq \\
        --in2 ${prefix}_2.fq \\
        --out1 ${prefix}_1.trim.fq \\
        --out2 ${prefix}_2.trim.fq \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --thread $task.cpus \\
        --detect_adapter_for_pe \\
        2> ${prefix}.fastp.log
    """
}

process INDEX {

    tag "$transcriptome"

    input:
    path transcriptome

    output:
    path 'index', emit: index


    script:
    """
    mkdir index
    salmon index --threads $task.cpus -t $transcriptome -i index
    """

}


process QUANT {

    tag "$sample_id"

    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    path index
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${prefix}"), emit: results
 
    script:
    prefix   = "${sample_id}"
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $prefix
    """
}


process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode:'copy'
       
    input:
    path (quant)
    path (fastqc)
    path (fastp)
    
    output:
    path 'multiqc_report.html'
     
    script:
    """
    multiqc . 
    """
}