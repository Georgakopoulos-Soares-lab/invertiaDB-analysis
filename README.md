# InvertiaDB Analysis & Extraction

## Introduction

MINDI stands for Mirror (MR), Inverted (IR), Direct (DI) repeats.

Mindi is a python package, and a wrapper around non b-gfa software,
that provides a python interface to extract non B-DNA motifs across any genome.

Mindi is using the snakemake workflow language (version 7.7.4).

We used the computational workflow to systematically extract, pre-process and validate inverted repeats for invertiaDB project.

In particular, in order to extract the inverted repeats, we used the following command:

```
./gfa -seq <infile> -out <outfile> \
                -minIR 8 \
                -maxSpacer 100 \
                -skipTriplex \
                -skipSlipped \
                -skipSTR \
                -skipWGET \
                -skipMR \
                -skipCruciform \
                -skipDR \
                -skipAPR \
                -skipZ \
                -skipGQ
```

The output file of this command is pre-processed in the following steps:

- The resulting command produces two files: a TSV and a GFF file.
- Each of these files contains the extractions of the inverted repeat sequences.
- Only the TSV file is being kept and is being used to extract the spacer sequence, the arm sequence 
and, finally, validate the sequence with the original fasta sequence of the genomic assembly.
- The file is being processed to contain only the inverted repeat sequences with arm length strictly more or equal to 10 base pairs and 
spacer length less or equal to 8 base pairs long.
- This results in a processed inverted repeat file, containing IR sequences with at least 10bp arms long and spacer at most 8bp long.

## Installation

### Local Installation

You can simply install it by cloning the directory locally:

```
git clone git@github.com:Georgakopoulos-Soares-lab/invertiaDB-analysis.git 
```

and then you can simply install:

```
pip install .
```

## Usage

mindi provides several pipelines for the extraction and analysis of non-b DNA motifs.

The workflows provided, are designed with big data in mind. The scheduler assigns a large chunk of genomes, 
provided by the design file (see below), into buckets. Subsequently, these buckets are being processed in parallel,
but each individual bucket and its associated genomic assemblies will be processed on a single core. 
The file assigned to buckets, takes place with a greedy strategy which optimizes the time efficiency on each individual core.

Notably, for the invertiaDB paper we extracted IR motifs across the complete genomes of
NCBI genbank & refseq databases. Note, that the input genomes don't necessarily have to be unzipped. 
The workflow is able to handle compressed files with the gzip strategy as well as raw fasta.

The pipelines are dependant on the scripts provided by mindi.

It also provides other utilities and pipelines.

Given a set of fasta files, we can extract the inverted repeats by providing:

- a design.csv file
- a configuration yaml file
- submitting our snakemake pipeline: `nonbdna_pipe.smk`

The design file holds the phylogenetic information of each assembly and is required step, 
before executing the pipeline. Note, that the phylogenetic information could be potentially empty,
but it is strongly encouraged to be present.

The configuration yaml file should look as follows:

```
buckets: 2
pattern: ['MR']
out: 'repeat_out'
suffix: 'fna'
DESIGN: 'design.csv'
step_threshold: 0.05
log_sleep: 0.5
```

The parameters are explained by the following table:

| Parameter | Explanation |
| ----------| ----------- |
| buckets   | total number of nodes used |
| pattern   | repeat pattern to extract, i.e. `IR`|
| out       | directory where pipeline stores the outputs |
| suffix    | suffix used to determine correct files |
| DESIGN    | file containing the accessions to extract |
| step_threshold | parameter used on GC threshold when estimating densities of repeats |
| log_sleep      | sleep parameter for daemon logging thread |

#### Core Dependencies

| Dependencies |
| -----------  |
| snakemake==7.4.3 |
| bedtools         |
| pybedtools       |

#### Extraction using snakemake

Finally, we need a directory of genomes:

```
├── GCA_000003135.1_ASM313v1_genomic.fna
├── GCA_000003215.1_ASM321v1_genomic.fna
├── GCA_000003645.1_ASM364v1_genomic.fna
├── GCA_000003925.1_ASM392v1_genomic.fna
├── GCA_000003955.1_ASM395v1_genomic.fna
```

To submit the pipeline, we simply execute the bash script:

```
export CORES=2
bash nonbdna_sub.sh $CORES
```

The command invokes the snakemake pipeline with the appropriate parameters as inherited from the configuration yaml files.

```
snakemake --snakefile nonbdna_pipe.smk \
            --configfile config/config.server.yaml \
	        --rerun-triggers mtime \
            --rerun-incomplete \
            --reason \
            --keep-going \
            --jobs $j \
            --latency-wait ${latency} \
            --cluster-config config/cluster_nonbdna.yaml \
            --cluster "sbatch -p {cluster.partition} \
                -t {cluster.time} \
                --mem={cluster.mem} \
                --nodes={cluster.nodes} \
                -J {cluster.jobName} \
                -o MindiNonBDNAJobDetails/{cluster.jobName}-%x-%j.out \
                -e MindiNonBDNAJobDetails/{cluster.jobName}-%x-%j.err"
```

Immediatly, the script understands if we are connected locally or not.
Note, that if you can to execute the script in a HPC server, an additional cluster.yaml 
is required, specifying the memory and time requirements of the HPC nodes.

## Docker
