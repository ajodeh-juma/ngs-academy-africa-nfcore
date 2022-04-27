#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// define workflow parameters
params.reads = "$projectDir/data/ggal/*_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.outdir = "results"




workflow {

    // data channel for transcriptome file
    ch_transcriptome = file(params.transcriptome_file, checkIfExists: true)
    
    // data channel for the input reads
    ch_raw_reads = Channel
        .fromFilePairs( params.reads, checkIfExists: true)
        .map { row ->
        def meta = [:]
        meta.id = row[0]
        [ meta , [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ]}


    // run fastqc process
    FASTQC(ch_raw_reads)
    ch_fastqc_multiqc = FASTQC.out.zip

    // run salmon processes
    INDEX(ch_transcriptome)
    QUANT(INDEX.out.index, ch_raw_reads)
    ch_salmon_multiqc = QUANT.out.results

    // run multiqc process
    MULTIQC( 
        ch_fastqc_multiqc.collect{it[1]}.ifEmpty([]),
        ch_salmon_multiqc.collect{it[1]}.ifEmpty([]),
    )
}

process FASTQC {

    tag "$meta.id"

    publishDir "${params.outdir}/fastqc", mode:'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip

    script:
    def prefix   = "${meta.id}"
    // Add soft-links to original FastQs for consistent naming in pipeline
    """
    [ ! -f  ${prefix}_1.fq ] && ln -s ${reads[0]} ${prefix}_1.fq
    [ ! -f  ${prefix}_2.fq ] && ln -s ${reads[1]} ${prefix}_2.fq

    fastqc --threads $task.cpus ${prefix}_1.fq ${prefix}_2.fq
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

    tag "$meta.id"

    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    path index
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}"), emit: results
 
    script:
    prefix   = "${meta.id}"
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $prefix
    """
}


process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode:'copy'
       
    input:
    path (quant)
    path (fastqc)
    
    output:
    path 'multiqc_report.html'
     
    script:
    """
    multiqc . 
    """
}


