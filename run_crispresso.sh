#!/bin/bash

annotation_file="/Users/smisek/OneDrive/Postdoc/POMT_and_Circumstance/sequencing_competition_assay/12.20.24_h9-NSC_E1/sample_annotation.txt"
number_of_primer_pairs=$(wc -l < "${annotation_file}" | tr -d '[:space:]')
input_dir="/Users/smisek/OneDrive/Postdoc/POMT_and_Circumstance/sequencing_competition_assay/12.20.24_h9-NSC_E1/t0/split_reads/"

# Validate input directory
if [ ! -d "${input_dir}" ]; then
    echo "Error: Input directory ${input_dir} does not exist."
    exit 1
fi

# Loop through all primer pairs and run CRISPResso on each file
for i in $(seq 1 $number_of_primer_pairs); do
    # Extract data from the annotation file
    line=$(sed -n "${i}p" "${annotation_file}")
    pair_name=$(echo "$line" | awk '{print $1}')
    amplicon_seq=$(echo "$line" | awk '{print $5}')
    guide_seq=$(echo "$line" | awk '{print $4}')

    # Validate variables
    if [ -z "$pair_name" ] || [ -z "$amplicon_seq" ] || [ -z "$guide_seq" ]; then
        echo "Error: Missing data for primer pair at line $i"
        continue
    fi

    #Now find all files that used that primer pair
    cd "${input_dir}"
    find ${pair_name}* > files_to_run.txt
    num_fastq=$(wc -l < files_to_run.txt | tr -d '[:space:]')

    if [ "${num_fastq}" -eq 0 ]; then
    echo "Error: No files found for primer pair ${pair_name}"
    continue
    fi

    for j in $(seq 1 ${num_fastq}); do 
        fastq=$(sed -n "${j}p" files_to_run.txt | awk '{print $1}')
        fastq_basename=$(basename "$fastq")
        output_folder="crispresso_output_${fastq_basename}"

        #Skip if not enough reads
        num_lines=$(wc -l < "${fastq}" | tr -d '[:space:]')
        if [ "${num_lines}" -lt 400 ]; then 
        echo "Warning: ${fastq} has fewer than 100 variants, skipping"
        continue
        fi

        docker run --platform linux/amd64 -v "${input_dir}:/DATA" \
            -w /DATA \
            -i pinellolab/crispresso2 CRISPResso \
            -g "${guide_seq}" \
            -r1 "${fastq_basename}" \
            --amplicon_seq "${amplicon_seq}" \
            -n "${output_folder}"
    done
done

