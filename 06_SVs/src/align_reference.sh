#!/usr/bin/env bash

# Info for S288C reference genome R65-3-1
ref_url='http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases'
ref_tgz='S288C_reference_genome_R64-3-1_20210421.tgz'
ref_fasta='S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa.gz'
ref_unzipped='S288C-reference.fasta'
sha256_desired='dbf065ffc3f5bbaa7554ef53c6861399eecc597fe058df5f81ba557944e3f86b'

# Download reference if it doesn't already exist
if [ ! -f "${ref_unzipped}" ]; then
    if [ ! -f "${ref_tgz}" ]; then
        wget -O ${ref_tgz} ${ref_url}/${ref_tgz}
        tar -zxvf  ${ref_tgz} --strip-components=1 ${ref_fasta} 
        gunzip -c $(basename ${ref_fasta}) > ${ref_unzipped}
    else
        echo "${ref_tgz} already downloaded"
    fi
fi

# Remove intermediates if they exist
rm -f -- ${ref_tgz}
rm -f -- $(basename ${ref_fasta})


# Ensure file is intact
sha256_actual=$(sha256sum ${ref_unzipped} | awk '{print $1}')
if [ "${sha256_desired}" != "${sha256_actual}" ]; then
    echo "checksum does not match, ${ref_unzipped} contents not as expected. Aborting."
    exit 0
fi

assemblies=($(ls assemblies/*.fasta))

module load parallel
module load minimap2/2.24
module load samtools/1.16.1

export out_dir='./reference_alignment'
mkdir -p ${out_dir}

# Align assemblies to reference
minimap_reference() {
    local ref_path=${1}
    local query_path=${2}

    local ref=$(basename ${ref_path%.fasta})
    local query=$(basename ${query_path%.fasta})

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

export -f minimap_reference

parallel -j 1 minimap_reference ::: ${ref_unzipped} ::: ${assemblies[@]}
