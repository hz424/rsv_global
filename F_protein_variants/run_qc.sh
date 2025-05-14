#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --job-name=flye_rsv_each
#SBATCH --time=3:00:00
#SBATCH --partition=devel
#SBATCH --account=kuin0129
#SBATCH --output=flye_rsv_each.%j.out
#SBATCH --error=flye_rsv_each.%j.err
#SBATCH --mem=60G

eval "$(conda shell.bash hook)"
conda activate ~/envs/rsv_env

THREADS=${SLURM_NTASKS}
GENOME_SIZE="15k"
OUT_BASE="flye_out"
FILTERED_DIR="filtered_reads"

mkdir -p "${OUT_BASE}" "${FILTERED_DIR}"

for FQ in reads/*.fastq; do
  SAMPLE=$(basename "${FQ}" .fastq)
  OUTDIR="${OUT_BASE}/${SAMPLE}"
  ASMF="${OUTDIR}/assembly.fasta"
  FILTERED_FQ="${FILTERED_DIR}/${SAMPLE}.fq"

  # skip if already assembled
  if [[ -f "${ASMF}" ]]; then
    echo "=== SKIP ${SAMPLE}: assembly already exists ==="
    continue
  fi

  echo "=== Processing ${SAMPLE} ==="

  # only filter if we haven't done it already
  if [[ -f "${FILTERED_FQ}" ]]; then
    echo "+++ FILTERED FASTQ exists, skipping NanoFilt for ${SAMPLE}"
  else
    echo ">>> Running NanoFilt for ${SAMPLE}"
    NanoFilt -q 7 -l 300 < "${FQ}" > "${FILTERED_FQ}"
  fi

done
