#!/bin/bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_crispresso.sh -a <annotation_file> -i <input_dir> [-d docker_image] [-p platform]

Arguments:
  -a  Path to the sample annotation file (required)
  -i  Directory containing split FASTQ files (required)
  -d  Docker image for CRISPResso (default: pinellolab/crispresso2)
  -p  Platform passed to docker run --platform (default: linux/amd64)
EOF
}

annotation_file=""
input_dir=""
docker_image="pinellolab/crispresso2"
docker_platform="linux/amd64"

while getopts ":a:i:d:p:h" opt; do
    case "${opt}" in
        a) annotation_file="${OPTARG}" ;;
        i) input_dir="${OPTARG}" ;;
        d) docker_image="${OPTARG}" ;;
        p) docker_platform="${OPTARG}" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "${annotation_file}" || -z "${input_dir}" ]]; then
    usage
    exit 1
fi

if [[ ! -f "${annotation_file}" ]]; then
    echo "Error: Annotation file ${annotation_file} does not exist."
    exit 1
fi

if [[ ! -d "${input_dir}" ]]; then
    echo "Error: Input directory ${input_dir} does not exist."
    exit 1
fi

input_dir="${input_dir%/}"
number_of_primer_pairs=$(wc -l < "${annotation_file}" | tr -d '[:space:]')

# Loop through all primer pairs and run CRISPResso on each file
for i in $(seq 1 "${number_of_primer_pairs}"); do
    line=$(sed -n "${i}p" "${annotation_file}")
    pair_name=$(echo "$line" | awk '{print $1}')
    amplicon_seq=$(echo "$line" | awk '{print $5}')
    guide_seq=$(echo "$line" | awk '{print $4}')

    if [[ -z "${pair_name}" || -z "${amplicon_seq}" || -z "${guide_seq}" ]]; then
        echo "Error: Missing data for primer pair at line ${i}"
        continue
    fi

    files_list=$(mktemp)

    find "${input_dir}" -maxdepth 1 -type f -name "${pair_name}*" -print > "${files_list}"
    num_fastq=$(wc -l < "${files_list}" | tr -d '[:space:]')

    if [[ "${num_fastq}" -eq 0 ]]; then
        echo "No files found for primer pair ${pair_name}. This is not necessarially a bad thing."
        rm -f "${files_list}"
        continue
    fi

    for j in $(seq 1 "${num_fastq}"); do
        fastq=$(sed -n "${j}p" "${files_list}" | awk '{print $1}')
        fastq_basename=$(basename "$fastq")
        clean_base="${fastq_basename%.fastq}"
        clean_base="${clean_base#CRISPResso_on_crisprseq_}"
        clean_base="${clean_base#crisprseq_}"
        output_folder="${clean_base}"

        num_lines=$(wc -l < "${fastq}" | tr -d '[:space:]')
        if [[ "${num_lines}" -lt 400 ]]; then
            echo "Warning: ${fastq} has fewer than 100 variants, skipping"
            continue
        fi

        docker run --platform "${docker_platform}" -v "${input_dir}:/DATA" \
            -w /DATA \
            -i "${docker_image}" CRISPResso \
            -g "${guide_seq}" \
            -r1 "${fastq_basename}" \
            --amplicon_seq "${amplicon_seq}" \
            -n "${output_folder}" \
            -wc 1 \
            -w 10 \
            --exclude_bp_from_left 0 \
            --exclude_bp_from_right 0

    done
    rm -f "${files_list}"
done
