#!/usr/bin/env bash
# conda activate /home/kunet.ae/ku1417/envs/rsv_env
module load parallel/20250422
set -euo pipefail

### CONFIG
REF_A="reference/RSV_A_F_cds.fasta"
REF_B="reference/RSV_B_F_cds.fasta"
REF_A_NAME=$(grep ">" "$REF_A" | head -n 1 | sed 's/>//' | awk '{print $1}')
REF_B_NAME=$(grep ">" "$REF_B" | head -n 1 | sed 's/>//' | awk '{print $1}')

SNPEFF_DB_A="RSV_A_F"
SNPEFF_DB_B="RSV_B_F"
THREADS=8                                 # threads used inside each job (8 threads for faster alignment)

# --- STRINGENT QC THRESHOLDS ---
DP_THRESH=20                              # Minimum depth for consensus calling (increased from 10)
MIN_MAPQ=20                               # Minimum mapping quality for reads
MIN_VARIANT_QUAL=100                      # Minimum variant quality (increased from 20)
MIN_VARIANT_DP=20                         # Minimum depth for variant calling (increased from 10)
MIN_COVERAGE_FRACTION=0.98                # Minimum coverage fraction for final analysis (increased from 0.95)
MAX_N_FRACTION=0.02                       # Maximum fraction of Ns in consensus (2%)
MIN_TYPING_RATIO=3.0                      # Minimum A:B or B:A ratio for confident typing

OUTDIR="workflow_run_AB"
mkdir -p "${OUTDIR}"/{bam,vcf,consensus,protein,ann,logs,tmp}

ASSIGNMENT_COVFILE="${OUTDIR}/assignment_coverage.tsv"
LOCKFILE="${OUTDIR}/.assign_cov.lock"
echo -e "Sample\tAssignedType\tMappedReads_A\tMappedReads_B\tReferenceLength\tMeanDepth\tCoveredBases\tCoverageFraction" > "$ASSIGNMENT_COVFILE"

echo "Indexing references..."
samtools faidx "$REF_A"
samtools faidx "$REF_B"
minimap2 -d "${REF_A}.mmi" "$REF_A"
minimap2 -d "${REF_B}.mmi" "$REF_B"
echo "Indexing complete."

process_sample() {
  local FQ="$1"
  local SAMPLE
  SAMPLE=$(basename "$FQ" .fq)
  local LOGFILE="${OUTDIR}/logs/${SAMPLE}.log"

  exec > "$LOGFILE" 2>&1
  echo "========================================"
  echo "=== Processing: $SAMPLE"
  echo "=== Date: $(date)"
  echo "========================================"

  local REF_A="$2" REF_B="$3" REF_A_NAME="$4" REF_B_NAME="$5"
  local SNPEFF_DB_A="$6" SNPEFF_DB_B="$7" THREADS="$8" OUTDIR="$9"
  local ASSIGNMENT_COVFILE="${10}"
  local TMPDIR="${OUTDIR}/tmp/${SAMPLE}"
  mkdir -p "$TMPDIR"

  local ASSIGNED_TYPE="" ASSIGNED_REF="" ASSIGNED_REF_NAME="" ASSIGNED_BAM="" ASSIGNED_SNPEFF_DB=""
  local FINAL_BAM_PATH="${OUTDIR}/bam/${SAMPLE}.assigned.bam"

  # --- 1) Alignment for Assignment ---
  echo "[$(date +%T)] → Aligning to Ref A: $REF_A_NAME"
  minimap2 -ax map-ont -t "$THREADS" "${REF_A}.mmi" "$FQ" \
    | samtools view -b -F 4 -@ "$THREADS" -o "${TMPDIR}/${SAMPLE}.A.mapped.bam" -
  echo "[$(date +%T)] → Aligning to Ref B: $REF_B_NAME"
  minimap2 -ax map-ont -t "$THREADS" "${REF_B}.mmi" "$FQ" \
    | samtools view -b -F 4 -@ "$THREADS" -o "${TMPDIR}/${SAMPLE}.B.mapped.bam" -

  local mapped_A mapped_B
  mapped_A=$(samtools view -c "${TMPDIR}/${SAMPLE}.A.mapped.bam")
  mapped_B=$(samtools view -c "${TMPDIR}/${SAMPLE}.B.mapped.bam")
  echo "[$(date +%T)]   Mapped read count A: $mapped_A"
  echo "[$(date +%T)]   Mapped read count B: $mapped_B"

  # --- 2) Assign Type (matches Methods text: higher read count; ties → Undetermined) ---
  if [[ "$mapped_A" -gt "$mapped_B" ]]; then
    ASSIGNED_TYPE="A"; ASSIGNED_REF="$REF_A"; ASSIGNED_REF_NAME="$REF_A_NAME"; ASSIGNED_SNPEFF_DB="$SNPEFF_DB_A"
    echo "[$(date +%T)] → Assigned Type: A; re-aligning to A for downstream"
    minimap2 -ax map-ont -t "$THREADS" "${REF_A}.mmi" "$FQ" \
      | samtools view -b -@ "$THREADS" - \
      | samtools sort -@ "$THREADS" -o "$FINAL_BAM_PATH" -
  elif [[ "$mapped_B" -gt "$mapped_A" ]]; then
    ASSIGNED_TYPE="B"; ASSIGNED_REF="$REF_B"; ASSIGNED_REF_NAME="$REF_B_NAME"; ASSIGNED_SNPEFF_DB="$SNPEFF_DB_B"
    echo "[$(date +%T)] → Assigned Type: B; re-aligning to B for downstream"
    minimap2 -ax map-ont -t "$THREADS" "${REF_B}.mmi" "$FQ" \
      | samtools view -b -@ "$THREADS" - \
      | samtools sort -@ "$THREADS" -o "$FINAL_BAM_PATH" -
  else
    ASSIGNED_TYPE="Undetermined"
    # Write to coverage file (flock removed - simple append is atomic for small writes)
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$SAMPLE" "$ASSIGNED_TYPE" "$mapped_A" "$mapped_B" "NA" "NA" "NA" "NA" >> "$ASSIGNMENT_COVFILE"
    rm -rf "$TMPDIR"
    echo "[$(date +%T)] === Finished: $SAMPLE (Undetermined) ==="
    return 0
  fi

  ASSIGNED_BAM="$FINAL_BAM_PATH"
  samtools index "$ASSIGNED_BAM"
  # remove assignment BAMs to save space
  rm -f "${TMPDIR}/${SAMPLE}.A.mapped.bam" "${TMPDIR}/${SAMPLE}.B.mapped.bam"
  rm -rf "$TMPDIR"

  # --- 3) Coverage (streamed; no depth cap) ---
  echo "[$(date +%T)] → Computing coverage vs $ASSIGNED_REF_NAME"
  local REF_LEN
  REF_LEN=$(awk '$0 !~ />/{L+=length($0)} END{print L+0}' "$ASSIGNED_REF")
  read mean cov <<<"$(
    samtools depth -aa -d 0 -r "$ASSIGNED_REF_NAME" "$ASSIGNED_BAM" \
    | awk -v len="$REF_LEN" '{sum+=$3; if($3>0) cov++} END{ if(len>0) printf("%.2f\t%d", sum/len, cov); else printf("0.00\t0") }'
  )"
  local frac
  frac=$(awk -v c="$cov" -v l="$REF_LEN" 'BEGIN{ if(l>0) printf("%.3f", c/l); else print "0.000" }')
  echo "[$(date +%T)]   MeanDepth=${mean}x | CoveredBases=${cov}/${REF_LEN} | Fraction=${frac}"

  # Write to coverage file (flock removed - simple append is atomic for small writes)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE" "$ASSIGNED_TYPE" "$mapped_A" "$mapped_B" "$REF_LEN" "$mean" "$cov" "$frac" >> "$ASSIGNMENT_COVFILE"

  # --- 4) Variant calling (stringent QC) ---
  echo "[$(date +%T)] → Calling variants (MAPQ≥${MIN_MAPQ})"
  bcftools mpileup -Ou --fasta-ref "$ASSIGNED_REF" --annotate FORMAT/DP --min-MQ "$MIN_MAPQ" --min-BQ 20 -d 0 "$ASSIGNED_BAM" \
    | bcftools call -m --ploidy 1 -v -Oz -o "${OUTDIR}/vcf/${SAMPLE}.assigned.raw.vcf.gz"
  bcftools index "${OUTDIR}/vcf/${SAMPLE}.assigned.raw.vcf.gz"

  # --- 5) Filter variants (stringent: QUAL≥100, DP≥20) ---
  echo "[$(date +%T)] → Filtering variants (QUAL≥${MIN_VARIANT_QUAL}, DP≥${MIN_VARIANT_DP})"
  bcftools filter --exclude "QUAL<${MIN_VARIANT_QUAL} || FORMAT/DP<${MIN_VARIANT_DP}" \
    "${OUTDIR}/vcf/${SAMPLE}.assigned.raw.vcf.gz" \
    -Oz -o "${OUTDIR}/vcf/${SAMPLE}.assigned.filtered.vcf.gz"
  bcftools index "${OUTDIR}/vcf/${SAMPLE}.assigned.filtered.vcf.gz"

  # --- 6) Build a BED mask of low-depth sites (DP < $DP_THRESH) and generate consensus with masking ---
  echo "[$(date +%T)] → Masking low-depth sites (<$DP_THRESH) to N and generating consensus (SNPs only)"
  local CONSENSUS_FASTA="${OUTDIR}/consensus/${SAMPLE}.consensus.fasta"
  local CHAIN_FILE="${OUTDIR}/consensus/${SAMPLE}.chain.txt"
  local FILTERED_VCF="${OUTDIR}/vcf/${SAMPLE}.assigned.filtered.vcf.gz"
  local FILTERED_SNP_VCF="${OUTDIR}/vcf/${SAMPLE}.assigned.filtered.snps.vcf.gz"
  local MASK="${OUTDIR}/consensus/${SAMPLE}.mask.bed"

  samtools depth -aa -d 0 -r "$ASSIGNED_REF_NAME" "$ASSIGNED_BAM" \
    | awk -v t="$DP_THRESH" '$3<t {printf "%s\t%d\t%d\n",$1,$2-1,$2}' > "$MASK"

  # Extract only SNPs (exclude indels to avoid frameshifts)
  bcftools view -v snps "$FILTERED_VCF" -Oz -o "$FILTERED_SNP_VCF"
  bcftools index "$FILTERED_SNP_VCF"

  bcftools consensus -f "$ASSIGNED_REF" -H A -m "$MASK" \
    "$FILTERED_SNP_VCF" -c "$CHAIN_FILE" -o "$CONSENSUS_FASTA" \
    || { echo "[$(date +%T)] ERROR: bcftools consensus failed for $SAMPLE"; exit 1; }

  if [[ -s "$CONSENSUS_FASTA" ]]; then
    sed -i "1 s|^>.*$|>${SAMPLE}_F_${ASSIGNED_TYPE}_consensus|" "$CONSENSUS_FASTA"
  else
    echo "[$(date +%T)] WARNING: Empty consensus for $SAMPLE"
  fi

  # --- 7) Translate CDS ---
  echo "[$(date +%T)] → Translating consensus CDS"
  local PROTEIN_FASTA="${OUTDIR}/protein/${SAMPLE}.F.protein.fa"
  if [[ -s "$CONSENSUS_FASTA" ]]; then
    transeq -sequence "$CONSENSUS_FASTA" -outseq "$PROTEIN_FASTA" -frame 1 -clean \
      || { echo "[$(date +%T)] ERROR: transeq failed for $SAMPLE"; exit 1; }
  else
    : > "$PROTEIN_FASTA"
  fi

  # --- 8) Annotation (SnpEff) + 9) Extract AA changes ---
  local ANNOTATED_VCF="${OUTDIR}/ann/${SAMPLE}.assigned.ann.vcf"
  local MUTATIONS_TXT="${OUTDIR}/ann/${SAMPLE}.assigned.mutations.txt"
  if [[ $(bcftools view -H "$FILTERED_VCF" | wc -l) -gt 0 ]]; then
    echo "[$(date +%T)] → Annotating with snpEff ($ASSIGNED_SNPEFF_DB)"
    snpEff -noStats -v "$ASSIGNED_SNPEFF_DB" "$FILTERED_VCF" > "$ANNOTATED_VCF" || : > "$ANNOTATED_VCF"
    if [[ -s "$ANNOTATED_VCF" ]]; then
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%ANN]\n' "$ANNOTATED_VCF" > "$MUTATIONS_TXT" || : > "$MUTATIONS_TXT"
    else
      : > "$MUTATIONS_TXT"
    fi
  else
    : > "$ANNOTATED_VCF"; : > "$MUTATIONS_TXT"
  fi

  echo "[$(date +%T)] === Finished: $SAMPLE (Assigned: $ASSIGNED_TYPE) ==="
  echo "========================================"
}

export -f process_sample
export REF_A REF_B REF_A_NAME REF_B_NAME SNPEFF_DB_A SNPEFF_DB_B THREADS OUTDIR ASSIGNMENT_COVFILE
export DP_THRESH MIN_MAPQ MIN_VARIANT_QUAL MIN_VARIANT_DP MIN_COVERAGE_FRACTION MAX_N_FRACTION MIN_TYPING_RATIO

TOTAL_CORES=$(nproc)
NUM_JOBS=5  # 5 parallel jobs (5 jobs × 8 threads = 40 CPUs)
if [[ $NUM_JOBS -lt 1 ]]; then NUM_JOBS=1; fi

echo "Starting parallel processing with $NUM_JOBS jobs (each using up to $THREADS threads internally)..."
echo "Input FASTQ files from: filtered_reads/"
echo "Output directory: $OUTDIR"
echo "Log files in: ${OUTDIR}/logs/"

find filtered_reads -name '*.fq' | sort \
  | parallel --eta -j "$NUM_JOBS" process_sample {} "$REF_A" "$REF_B" "$REF_A_NAME" "$REF_B_NAME" "$SNPEFF_DB_A" "$SNPEFF_DB_B" "$THREADS" "$OUTDIR" "$ASSIGNMENT_COVFILE"

echo "[$(date +%T)] All processing complete."
echo "Assignment and coverage metrics saved to $ASSIGNMENT_COVFILE"

# --- 10) Mutation matrix (corrected - uses consensus sequences) ---
MATRIX="${OUTDIR}/F_mutation_matrix.tsv"
echo "[$(date +%T)] Building corrected mutation matrix from consensus sequences..."
python3 << 'EOF'
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import pandas as pd
import sys

# Load reference sequences
REF_A_CDS = "reference/RSV_A_F_cds.fasta"
REF_B_CDS = "reference/RSV_B_F_cds.fasta"
OUTDIR = "workflow_run_AB"
ASSIGNMENT_FILE = f"{OUTDIR}/assignment_coverage.tsv"
CONSENSUS_DIR = f"{OUTDIR}/consensus"
MATRIX_FILE = f"{OUTDIR}/F_mutation_matrix.tsv"
MIN_COVERAGE_FRACTION = 0.98  # Stringent: 98% coverage required
MAX_N_FRACTION = 0.02         # Stringent: max 2% Ns allowed

# Load reference proteins
ref_seqs = {}
for ref_file, ref_type in [(REF_A_CDS, 'A'), (REF_B_CDS, 'B')]:
    for record in SeqIO.parse(ref_file, "fasta"):
        ref_protein = str(record.seq.translate())
        ref_seqs[ref_type] = ref_protein

# Load assignment data
cov_df = pd.read_csv(ASSIGNMENT_FILE, sep='\t')
cov_df['CoverageFraction'] = pd.to_numeric(cov_df['CoverageFraction'], errors='coerce')
cov_df['MeanDepth'] = pd.to_numeric(cov_df['MeanDepth'], errors='coerce')

# Stringent QC filters
good_df = cov_df[
    (cov_df['CoverageFraction'] >= MIN_COVERAGE_FRACTION) &
    (cov_df['AssignedType'].isin(['A','B'])) &
    (cov_df['MeanDepth'] >= 20)  # Minimum 20x mean depth
]
good_samples = good_df.set_index('Sample')['AssignedType'].to_dict()

print(f"Samples passing QC (Cov≥{MIN_COVERAGE_FRACTION}, Depth≥20x): {len(good_samples)}", file=sys.stderr)

# Process each sample
sample_muts = {}
skipped_n_content = 0
skipped_premature_stop = 0

for sample, assigned_type in good_samples.items():
    consensus_file = Path(CONSENSUS_DIR) / f"{sample}.consensus.fasta"
    if not consensus_file.exists():
        continue

    try:
        record = next(SeqIO.parse(consensus_file, "fasta"))
        consensus_seq = str(record.seq)

        # Skip if too many Ns (stringent: 2%)
        n_count = consensus_seq.count('N') + consensus_seq.count('n')
        if n_count / len(consensus_seq) > MAX_N_FRACTION:
            skipped_n_content += 1
            continue

        # Translate consensus
        consensus_protein = str(Seq(consensus_seq).translate())
        ref_protein = ref_seqs[assigned_type]

        # Check for premature stop codons (excluding natural stop at end)
        premature_stop_pos = None
        for pos in range(len(consensus_protein) - 1):
            if consensus_protein[pos] == '*':
                premature_stop_pos = pos
                skipped_premature_stop += 1
                break

        # Skip samples with premature stops (indicates frameshift/poor quality)
        if premature_stop_pos is not None:
            continue

        # Find mutations
        mutations = set()
        max_pos = min(len(ref_protein), len(consensus_protein))

        for pos in range(max_pos):
            ref_aa = ref_protein[pos]
            cons_aa = consensus_protein[pos]
            if cons_aa == 'X' or ref_aa == cons_aa:
                continue
            mutations.add(f"{assigned_type}:{ref_aa}{pos+1}{cons_aa}")

        if mutations:
            sample_muts[sample] = mutations
    except:
        continue

print(f"Skipped {skipped_n_content} samples due to high N content (>{MAX_N_FRACTION*100}%)", file=sys.stderr)
print(f"Skipped {skipped_premature_stop} samples due to premature stop codons", file=sys.stderr)

# Build mutation matrix
all_muts = sorted({m for s in sample_muts.values() for m in s},
                  key=lambda x: (x.split(':')[0], int(''.join(filter(str.isdigit, x))), x))
cols = sorted(good_samples.keys())
df = pd.DataFrame(0, index=all_muts, columns=cols)
for sample, muts in sample_muts.items():
    for mut in muts:
        if mut in df.index:
            df.at[mut, sample] = 1
df.to_csv(MATRIX_FILE, sep='\t', index_label='Mutation')
print(f"Mutation matrix: {len(all_muts)} mutations, {len(sample_muts)} samples", file=sys.stderr)
EOF

echo "[$(date +%T)] → Mutation matrix generation complete."

# Remove old incorrect mutation matrix if it exists
OLD_MATRIX="${OUTDIR}/F_mutation_matrix_corrected.tsv"
if [[ -f "$OLD_MATRIX" ]]; then
  rm "$OLD_MATRIX"
  echo "[$(date +%T)] → Removed old temporary mutation matrix file"
fi

echo "[$(date +%T)] Workflow finished."
