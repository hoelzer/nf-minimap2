# Nextflow example using `minimap2`

A small example workflow to introduce Nextflow and Conda or Docker integration. 

First, install [Nextflow](https://nextflow.io/), [Conda](https://docs.conda.io/en/latest/miniconda.html), and [Docker](https://docs.docker.com/engine/installation/).

Second, get example data via:
```bash
wget https://raw.githubusercontent.com/KleistLab/GInPipe/main/demo/demo_reference.fasta .
wget https://raw.githubusercontent.com/KleistLab/GInPipe/main/demo/demo_samples.fasta .
```
(obtained from https://github.com/KleistLab/GInPipe/tree/main/demo, can be also found in the `data` folder in this repository)

Then run the workflow:
```bash
nextflow run main.nf --reference demo_reference.fasta --query demo_samples.fasta
```

The workflow will align the FASTA sequences in the query file vs. the target sequence in the reference FASTA file.

**Hint**: The workflow is configured to use some default input data (see the `nextflow.config` file). Therefore, you can also just run it via
```bash
nextflow run main.nf
```
