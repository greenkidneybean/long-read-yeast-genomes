#!/usr/bin/env bash

export ZIP_FILENAME=${1}
export PCT_IDENTITY=${2}
export NTHREADS=${3}
export E_VALUE_CUTOFF=${4}
export BLAST_DB=${5}

readarray -t FILES < <(zipinfo -1 ${ZIP_FILENAME})

echo ${FILES[@]}

function check_files() {
    local filename=${1}
    local cross=$(echo ${filename%.sam} | cut -d '_' -f 1)
    local plate=$(echo ${filename%.sam} | cut -d '_' -f 2)
    local well=$(echo ${filename%.sam} | cut -d '_' -f 3)

    echo $cross $plate $well

    local blast_file="data/blasts/${cross}_${plate}_${well}.blast.gz"
    local fasta_file="data/fasta/${cross}_${plate}_${well}.fasta"
    if [ ! -f "${blast_file}" ]; then
        if [ ! -f "${fasta_file}" ]; then
            echo "building  ${fasta_file}"
            build_fasta ${cross} ${plate} ${well}
        fi
        echo "building ${blast_file}"
        run_blast ${fasta_file} ${blast_file}
    fi
    if [ -f "${fasta_file}" ]; then
        echo "removing ${fasta_file}"
        rm ${fasta_file}
    fi
}

function build_fasta() {
    local cross=${1}
    local plate=${2}
    local well=${3}
    samtools fasta \
    -1 data/fasta/${cross}_${plate}_${well}.fasta \
    -2 data/fasta/${cross}_${plate}_${well}.fasta \
    -0 data/fasta/${cross}_${plate}_${well}.fasta \
    -s data/fasta/${cross}_${plate}_${well}.fasta \
    <(cat data/header-${cross}.txt <(unzip -p  ${ZIP_FILENAME} ${cross}_${plate}_${well}.sam))
}

function run_blast() {
    local fasta_file=${1}
    local blast_file=${2}
    ./blast/bin/blastn \
        -perc_identity ${PCT_IDENTITY} \
        -num_threads 1 \
        -evalue ${E_VALUE_CUTOFF} \
        -outfmt 6 \
        -query ${fasta_file} \
        -db ${BLAST_DB}  | gzip -c > ${blast_file}
}

export -f check_files
export -f build_fasta
export -f run_blast

parallel -j ${NTHREADS} check_files ::: ${FILES[@]}
