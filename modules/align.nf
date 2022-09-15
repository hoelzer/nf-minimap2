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
