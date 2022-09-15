#!/usr/bin/env nextflow

nextflow.enable.dsl=2

if (params.reference) { input_reference_fasta = Channel.fromPath(params.reference) } 
if (params.query) { input_query_fasta = Channel.fromPath(params.query) }

process ALIGN {

    // pull an image from Dockerhub or us already available local version
    container 'mhoelzer/minimap2:2.24'

    publishDir "results", mode: 'copy', pattern: "*.sam"

    input: 
    tuple path(query), path(reference)

    output:
    path("${query.simpleName}.sam")

    script:
    """
    minimap2 -ax asm5 ${reference} ${query} > ${query.simpleName}.sam
    """
}

workflow {
    ALIGN(input_query_fasta.combine(input_reference_fasta))
}