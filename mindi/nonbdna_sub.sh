#!/bin/bash

j=${1:-1}
latency=${2:-20}

if [[ ! -n "$SSH_CONNECTION" ]];
then
  echo "Local environment detected."
  echo "Initializing bioinformatics genomic compartment enrichment analysis. (mode ${mode}; cores ${j}) [LOCAL]. Authored by Nikol Chantzi <3."
	snakemake --snakefile nonbdna_pipe.smk \
            --configfile config/config.yaml \
	          --rerun-triggers mtime \
            --rerun-incomplete \
            --reason \
            --keep-going \
            --latency-wait ${latency} \
            --cores $j
else
  echo "SSH Connection detected."
  echo "Initializing bioinformatics genomic compartment enrichment analysis. (mode ${mode}; cores ${j}) [SERVER]. Authored by Nikol Chantzi <3."
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
fi
