#!/usr/bin/env bash

cd reference_alignment

for file in *.coords; do
    if head -n 1 ${file}  | grep -q 'LEN R'; then
        echo "${file} already fixed, skipping"
    else
        echo "trimming extra lines from top of ${file}"
        tail -n +4 ${file} > ${file}.tmp
        mv ${file}.tmp ${file}
    fi
done
