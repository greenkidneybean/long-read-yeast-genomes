#!/usr/bin/env bash

snakemake -k -j 20 --max-jobs-per-second 1 --max-status-checks-per-second 1 --use-singularity --latency-wait=180 --cluster="sbatch -c {threads} --mem={resources.mem_mb} --time={resources.runtime_min}"
