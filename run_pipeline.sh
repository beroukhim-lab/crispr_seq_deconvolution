#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'EOF'
Usage: run_pipeline.sh -f <fastq> | -F <fastq_dir> -a <annotation_file> -o <output_dir> [-d docker_image] [-p platform]

Runs:
  1) split_fastq.py   (-i FASTQ -a ANNOTATION -o OUTPUT_DIR)
  2) qual_check.py    (-i OUTPUT_DIR -a ANNOTATION)
  3) run_crispresso.sh (-a ANNOTATION -i OUTPUT_DIR [-d docker_image] [-p platform])
  4) summarize_crispresso_indels.py (-i OUTPUT_DIR_ROOT -o <fastq_root>/allele_deletion_summary.tsv)
EOF
}

fastq=""
fastq_dir=""
annotation=""
output_dir=""
docker_image="pinellolab/crispresso2"
docker_platform="linux/amd64"

while getopts ":f:F:a:o:d:p:h" opt; do
    case "${opt}" in
        f) fastq="${OPTARG}" ;;
        F) fastq_dir="${OPTARG}" ;;
        a) annotation="${OPTARG}" ;;
        o) output_dir="${OPTARG}" ;;
        d) docker_image="${OPTARG}" ;;
        p) docker_platform="${OPTARG}" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -n "${fastq}" && -n "${fastq_dir}" ]]; then
    echo "Error: specify either a single FASTQ (-f) or a directory of FASTQs (-F), not both."
    usage
    exit 1
fi

if [[ -z "${fastq}" && -z "${fastq_dir}" ]]; then
    echo "Error: must provide a FASTQ (-f) or a directory of FASTQs (-F)."
    usage
    exit 1
fi

if [[ -z "${annotation}" || -z "${output_dir}" ]]; then
    usage
    exit 1
fi
output_dir="${output_dir%/}"
fastq_root=""

if [[ ! -f "${annotation}" ]]; then
    echo "Error: Annotation file ${annotation} not found."
    exit 1
fi

fastq_list=()

if [[ -n "${fastq}" ]]; then
    if [[ ! -f "${fastq}" ]]; then
        echo "Error: FASTQ file ${fastq} not found."
        exit 1
    fi
    fastq_root="$(cd "$(dirname "${fastq}")" && pwd)"
    fastq_list+=("${fastq}")
else
    if [[ ! -d "${fastq_dir}" ]]; then
        echo "Error: FASTQ directory ${fastq_dir} not found."
        exit 1
    fi
    fastq_root="$(cd "${fastq_dir}" && pwd)"
    while IFS= read -r fq; do
        fastq_list+=("$fq")
    done < <(find "${fastq_dir}" -maxdepth 1 -type f -name "*.fastq" -print | sort)

    if [[ ${#fastq_list[@]} -eq 0 ]]; then
        echo "Error: no FASTQ files (*.fastq) found in ${fastq_dir}"
        exit 1
    fi
fi

for fq in "${fastq_list[@]}"; do
    base=$(basename "${fq}")
    base="${base%.fastq}"
    run_outdir="${output_dir}/${base}"
    mkdir -p "${run_outdir}"

    echo "=== Processing ${fq} ==="
    echo "Splitting FASTQ..."
    python3 "${SCRIPT_DIR}/split_fastq.py" -i "${fq}" -a "${annotation}" -o "${run_outdir}"

    echo "Running quality checks..."
    python3 "${SCRIPT_DIR}/qual_check.py" -i "${run_outdir}" -a "${annotation}"

    echo "Running CRISPResso..."
    bash "${SCRIPT_DIR}/run_crispresso.sh" -a "${annotation}" -i "${run_outdir}" -d "${docker_image}" -p "${docker_platform}"
done

summary_out="${fastq_root}/allele_deletion_summary.tsv"
echo "Summarizing CRISPResso indels into ${summary_out}..."
if ! python3 "${SCRIPT_DIR}/summarize_crispresso_indels.py" -i "${output_dir}" -o "${summary_out}"; then
    echo "Warning: could not write summary to ${summary_out}. Writing inside output_dir instead."
    summary_out="${output_dir}/allele_deletion_summary.tsv"
    python3 "${SCRIPT_DIR}/summarize_crispresso_indels.py" -i "${output_dir}" -o "${summary_out}"
fi

echo "Pipeline complete. Outputs written under ${output_dir}"
echo "Summary table: ${summary_out}"
echo "Stacked bar plot: ${summary_out%.tsv}_stacked_bar.pdf"
echo "Pie charts PDF: ${summary_out%.tsv}_pies.pdf"
