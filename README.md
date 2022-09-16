# Nextflow example using `minimap2`

A small example workflow to introduce `Nextflow` and `Conda` or `Docker` integration. We will use the famous mapping tool [`minimap2`](https://github.com/lh3/minimap2) as example tool. 

## Setup

First, install [Nextflow](https://nextflow.io/), [Conda](https://docs.conda.io/en/latest/miniconda.html), and [Docker](https://docs.docker.com/engine/installation/).

Second, clone this repository and change into the new dir:
```bash
git clone https://github.com/hoelzer/nf-minimap2.git
cd nf-minimap2
```

Third, get some example data via:
```bash
wget https://raw.githubusercontent.com/KleistLab/GInPipe/main/demo/demo_reference.fasta .
wget https://raw.githubusercontent.com/KleistLab/GInPipe/main/demo/demo_samples.fasta .
```
(obtained from https://github.com/KleistLab/GInPipe/tree/main/demo, can be also found in the `data` folder in this repository)

## Run

We will show examples of using `minimap2` directly on the command line in a `Conda` environment and within a small `Nextflow` pipeline. Within `Nextflow`, we will per default use a `Docker` container that holds `minimap2` and switch as an alternative also to `Conda`. 

### Manuallly

```bash
# create a new Conda environment and install minimap2
conda create -n minimap2 -c bioconda minimap2:2.24
# activate the new env
conda activate minimap2
# does the tool work?
minimap2 --help
minimap2 --version

# now we align some FASTA sequences against another reference sequence
# we use minimap2's option 'asm5' which is actually intended for cross-species full-genome alignment
# Attention: you likely need to switch to other parameters for Illumina or Nanopore reads! 
REFERENCE=demo_reference.fasta
QUERY=demo_samples.fasta
minimap2 -ax asm5 $REFERENCE $QUERY > $(basename $QUERY .fasta).sam

# as you can see, we define two variables to point to the two input files.
# we also use another command 'basename' to get the basename of the $QUERY variable
# we pipe the output into a new file with the file ending SAM
# see https://samtools.github.io/hts-specs/SAMv1.pdf for SAM specifications
```

When the calculation is done you can inspect the `SAM` file. What does it tell you? Did the alignment work? 

### Nextflow

Then run the workflow:
```bash
nextflow run main.nf --reference demo_reference.fasta --query demo_samples.fasta
```

The workflow will align the FASTA sequences in the query file vs. the target sequence in the reference FASTA file.

**Hint**: The workflow is configured to use some default input data (see the `nextflow.config` file). Therefore, you can also just run it via
```bash
nextflow run main.nf
```

#### Use Conda instead of Docker

Per default, the workflow runs with `Docker` support using the image defined in the `main.nf` file. However, you can also use `Conda`. To do so, go in the `main.nf` file and comment the container definition and un-comment the conda defintion pointing to the `envs/minimap2.yaml` file like that:

```java
process ALIGN {

    // pull an image from Dockerhub or us already available local version.
    //container 'mhoelzer/minimap2:2.24'

    // use a conda env. If it does not exist, it's autonatically created.
    conda 'envs/minimap2.yaml'
...
}
```

When you start the workflow now, `Nextflow` will check for an available `Conda` environment and if it does not exist, create it for you. It just need to be created once. 

#### Run the process in a separate module

Per default, the workflow integrates the process `ALIGN` in the `main.nf` but that is **not** best practice. It is better to modularize processes. See the `align.nf` file in the `modules` folder. You can switch to using that process by simply commenting the whole `ALIGN` process code block in the `main.nf` file and un-commenting the `include` statement like that:

```java
// [1] here we directly include the process in the main workflow file
/*
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
*/
// [2] but it's better to separate main workflow and processes:
include { ALIGN } from './modules/align.nf'
```

Then, run the workflow again. 
