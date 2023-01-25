#!/usr/bin/env bash

assemblies=($(ls assemblies/*.fasta))

module load parallel
module load minimap2/2.24
module load samtools/1.16.1
export out_dir='./pairwise_alignments'

mkdir -p ${out_dir}

minimap_pairwise() {
    local ref_path=${1}
    local query_path=${2}

    local ref=$(basename ${ref_path%.fasta})
    local query=$(basename ${query_path%.fasta})

    if [ "${ref}" == "${query}" ]; then
        echo "${ref} and ${query} are same, skipping" && return
    fi

    local output="${out_dir}/${ref}_${query}.paf"

    if [ ! -f "${output}" ]; then
        echo "Generating ${output}"
        minimap2 -csd --secondary no \
            ${ref_path} \
            ${query_path} > ${output}
        return
    else
        echo "${output} already exists, skipping" && return
    fi
}

export -f minimap_pairwise

parallel -j 1 minimap_pairwise ::: ${assemblies[@]} ::: ${assemblies[@]}
