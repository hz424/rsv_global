#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --job-name=filter_rsv_reads_each # Renamed job to reflect new purpose
#SBATCH --time=3:00:00
#SBATCH --partition=devel
#SBATCH --account=kuin0129
#SBATCH --output=filter_rsv_reads_each.%j.out # Adjusted output file name
#SBATCH --error=filter_rsv_reads_each.%j.err  # Adjusted error file name
#SBATCH --mem=60G

eval "$(conda shell.bash hook)"
conda activate ~/envs/rsv_env

THREADS=${SLURM_NTASKS}
# GENOME_SIZE="15k" # GENOME_SIZE was likely for Flye, can be removed if not used elsewhere
FILTERED_DIR="filtered_reads"

mkdir -p "${FILTERED_DIR}" # Only create the filtered reads directory

for FQ in reads/*.fastq; do
  SAMPLE=$(basename "${FQ}" .fastq)
  FILTERED_FQ="${FILTERED_DIR}/${SAMPLE}.fq"

  echo "=== Processing ${SAMPLE} for filtering ===" # Clarified processing step

  # only filter if we haven't done it already
  if [[ -f "${FILTERED_FQ}" ]]; then
    echo "+++ FILTERED FASTQ exists, skipping NanoFilt for ${SAMPLE}"
  else
    echo ">>> Running NanoFilt for ${SAMPLE}"
    NanoFilt -q 7 -l 300 < "${FQ}" > "${FILTERED_FQ}"
  fi

done