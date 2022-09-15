#!/usr/bin/env nextflow

nextflow.enable.dsl=2

if (params.reference) { input_reference_fasta = Channel.fromPath(params.reference) } 
if (params.query) { input_query_fasta = Channel.fromPath(params.query) }

// [1] here we directly include the process in the main workflow file
process ALIGN {

    // pull an image from Dockerhub or us already available local version.
    container 'mhoelzer/minimap2:2.24'

    // use a conda env. If it does not exist, it's autonatically created.
    //conda 'envs/minimap2.yaml'

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
// [2] but it's better to separate main workflow and processes:
//include { ALIGN } from './modules/align.nf'

workflow {
    ALIGN(input_query_fasta.combine(input_reference_fasta))
}