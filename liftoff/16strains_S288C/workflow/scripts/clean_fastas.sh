

# Extract all fasta into current folder
for file in *.fasta; do
    strain=$(echo $file | cut -d "_" -f 1)
    sed -i 's/with.*//g' ${file}
    sed -i "s/tig.*chr/${strain}_chr/g" ${file}
    sed -i "s/tig.*mito/${strain}_mito/g" ${file}
    sed -i 's/ .*$//g' ${file}
done
# rename files
for file in *.fasta; do
    strain=$(echo $file | cut -d "_" -f 1)
    mv ${file} ${strain}.fasta
done
# archive
ls *.fasta | tar -czv -f pacbio_renamed_fastas.tar.gz --files-from -
