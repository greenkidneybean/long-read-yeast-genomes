#!/usr/bin/env bash

export ZIP_FILENAME=${1}
export PCT_IDENTITY=${2}
export NTHREADS=${3}
export E_VALUE_CUTOFF=${4}
export BLAST_DB=${5}

readarray -t FILES < <(zipinfo -1 ${ZIP_FILENAME} | head -n 10)

echo ${FILES[@]}

function check_files() {
    local filename=${1}
    local cross=$(echo ${filename%.sam} | cut -d '_' -f 1)
    local plate=$(echo ${filename%.sam} | cut -d '_' -f 2)
    local well=$(echo ${filename%.sam} | cut -d '_' -f 3)

    echo $cross $plate $well

    #local blast_file="data/blasts/${cross}_${plate}_${well}.blast.gz"
    local fastq_file="data/fastq/${cross}_${plate}_${well}.fastq"
    if [ ! -f "${blast_file}" ]; then
        if [ ! -f "${fastq_file}" ]; then
            echo "building  ${fastq_file}"
            build_fastq ${cross} ${plate} ${well}
        fi
        # echo "building ${blast_file}"
        # run_blast ${fastq_file} ${blast_file}
    fi
    # if [ -f "${fastq_file}" ]; then
    #     echo "removing ${fastq_file}"
    #     rm ${fastq_file}
    # fi
}

function build_fastq() {
    local cross=${1}
    local plate=${2}
    local well=${3}
    samtools fastq \
    -1 data/fastq/${cross}_${plate}_${well}.fastq \
    -2 data/fastq/${cross}_${plate}_${well}.fastq \
    -0 data/fastq/${cross}_${plate}_${well}.fastq \
    -s data/fastq/${cross}_${plate}_${well}.fastq \
    <(cat data/header-${cross}.txt <(unzip -p  ${ZIP_FILENAME} ${cross}_${plate}_${well}.sam))
}

function run_blast() {
    local fastq_file=${1}
    local blast_file=${2}
    ./blast/bin/blastn \
        -perc_identity ${PCT_IDENTITY} \
        -num_threads 1 \
        -evalue ${E_VALUE_CUTOFF} \
        -outfmt 6 \
        -query ${fastq_file} \
        -db ${BLAST_DB}  | gzip -c > ${blast_file}
}

export -f check_files
export -f build_fastq
export -f run_blast

parallel -j ${NTHREADS} check_files ::: ${FILES[@]}
