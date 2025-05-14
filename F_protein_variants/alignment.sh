#!/usr/bin/env bash
# conda activate /home/kunet.ae/ku1417/envs/rsv_env
module load parallel/20250422
# Ensure necessary tools like samtools, bcftools, minimap2, emboss, snpEff are loaded or in PATH

set -euo pipefail

### CONFIG
REF_A="reference/RSV_A_F_cds.fasta"
REF_B="reference/RSV_B_F_cds.fasta"
# extract sequence names to avoid manual errors
REF_A_NAME=$(grep ">" "$REF_A" | head -n 1 | sed 's/>//' | awk '{print $1}')
REF_B_NAME=$(grep ">" "$REF_B" | head -n 1 | sed 's/>//' | awk '{print $1}')

# SnpEff databases
SNPEFF_DB_A="RSV_A_F"
SNPEFF_DB_B="RSV_B_F"
THREADS=20 # Threads per sample processing job

OUTDIR="workflow_run_AB"
mkdir -p "${OUTDIR}"/{bam,vcf,consensus,protein,ann,logs,tmp}

# Combined Coverage/Assignment metrics file
ASSIGNMENT_COVFILE="${OUTDIR}/assignment_coverage.tsv"
echo -e "Sample\tAssignedType\tMappedReads_A\tMappedReads_B\tReferenceLength\tMeanDepth\tCoveredBases\tCoverageFraction" > "$ASSIGNMENT_COVFILE"

# --- Pre-index references ---
echo "Indexing references..."
samtools faidx "$REF_A"
samtools faidx "$REF_B"
minimap2 -d "${REF_A}.mmi" "$REF_A"
minimap2 -d "${REF_B}.mmi" "$REF_B"
echo "Indexing complete."

# --- Function to process a single sample ---
process_sample() {
  # Define local variables for clarity and safety
  local FQ="$1"
  local SAMPLE
  SAMPLE=$(basename "$FQ" .fq) # Extract sample name from fastq filename
  local LOGFILE="${OUTDIR}/logs/${SAMPLE}.log"

  # Redirect all output (stdout & stderr) for this function run to the log file
  exec > "$LOGFILE" 2>&1

  echo "========================================"
  echo "=== Processing: $SAMPLE"
  echo "=== Date: $(date)"
  echo "========================================"

  # Pass global variables explicitly to the function scope
  local REF_A="$2"
  local REF_B="$3"
  local REF_A_NAME="$4"
  local REF_B_NAME="$5"
  local SNPEFF_DB_A="$6"
  local SNPEFF_DB_B="$7"
  local THREADS="$8"
  local OUTDIR="$9"
  local ASSIGNMENT_COVFILE="${10}"
  local TMPDIR="${OUTDIR}/tmp/${SAMPLE}" # Temporary directory for intermediate files
  mkdir -p "$TMPDIR"

  # Variables to store assignment results
  local ASSIGNED_TYPE=""
  local ASSIGNED_REF=""
  local ASSIGNED_REF_NAME=""
  local ASSIGNED_BAM=""
  local ASSIGNED_SNPEFF_DB=""
  local FINAL_BAM_PATH="${OUTDIR}/bam/${SAMPLE}.assigned.bam" # Final sorted BAM path

  # --- 1) Alignment for Assignment ---
  echo "[$(date +%T)] → Aligning to Ref A: $REF_A_NAME"
  minimap2 -ax map-ont -t "$THREADS" "${REF_A}.mmi" "$FQ" \
    | samtools view -b -F 4 -@ "$THREADS" -o "${TMPDIR}/${SAMPLE}.A.mapped.bam" -

  echo "[$(date +%T)] → Aligning to Ref B: $REF_B_NAME"
  minimap2 -ax map-ont -t "$THREADS" "${REF_B}.mmi" "$FQ" \
    | samtools view -b -F 4 -@ "$THREADS" -o "${TMPDIR}/${SAMPLE}.B.mapped.bam" -

  local mapped_A=$(samtools view -c "${TMPDIR}/${SAMPLE}.A.mapped.bam")
  local mapped_B=$(samtools view -c "${TMPDIR}/${SAMPLE}.B.mapped.bam")
  echo "[$(date +%T)]   Mapped reads count A: $mapped_A"
  echo "[$(date +%T)]   Mapped reads count B: $mapped_B"

  # --- 2) Assign Type ---
  if [[ "$mapped_A" -gt "$mapped_B" ]]; then
    ASSIGNED_TYPE="A"
    ASSIGNED_REF="$REF_A"
    ASSIGNED_REF_NAME="$REF_A_NAME"
    ASSIGNED_SNPEFF_DB="$SNPEFF_DB_A"
    echo "[$(date +%T)] → Assigned Type: A (using A alignment for downstream analysis)"
    minimap2 -ax map-ont -t "$THREADS" "${REF_A}.mmi" "$FQ" \
        | samtools view -b -@ "$THREADS" - \
        | samtools sort -@ "$THREADS" -o "$FINAL_BAM_PATH" -
    ASSIGNED_BAM="$FINAL_BAM_PATH"

  elif [[ "$mapped_B" -gt "$mapped_A" ]]; then
    ASSIGNED_TYPE="B"
    ASSIGNED_REF="$REF_B"
    ASSIGNED_REF_NAME="$REF_B_NAME"
    ASSIGNED_SNPEFF_DB="$SNPEFF_DB_B"
    echo "[$(date +%T)] → Assigned Type: B (using B alignment for downstream analysis)"
    minimap2 -ax map-ont -t "$THREADS" "${REF_B}.mmi" "$FQ" \
        | samtools view -b -@ "$THREADS" - \
        | samtools sort -@ "$THREADS" -o "$FINAL_BAM_PATH" -
    ASSIGNED_BAM="$FINAL_BAM_PATH"
  else
    ASSIGNED_TYPE="Undetermined"
    echo "[$(date +%T)] → WARNING: Could not determine type (A=$mapped_A, B=$mapped_B). Skipping downstream analysis for this sample."
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$SAMPLE" "$ASSIGNED_TYPE" "$mapped_A" "$mapped_B" "NA" "NA" "NA" "NA" >> "$ASSIGNMENT_COVFILE"
    rm -rf "$TMPDIR"
    echo "[$(date +%T)] === Finished: $SAMPLE (Undetermined) ==="
    return 0
  fi

  samtools index "$ASSIGNED_BAM"
  rm -rf "$TMPDIR"

  # --- 3) Compute Coverage (on assigned BAM) ---
  echo "[$(date +%T)] → Computing coverage against assigned reference: $ASSIGNED_REF_NAME"
  local REF_LEN=$(awk '$0 !~ />/{total_len += length($0)} END{print total_len}' "$ASSIGNED_REF")
  local cov_data
  cov_data=$(samtools depth -a -r "$ASSIGNED_REF_NAME" "$ASSIGNED_BAM")

  local mean=0
  local cov=0
  if [[ -n "$cov_data" ]]; then
      read mean cov <<<$(printf '%s\n' "$cov_data" \
        | awk -v len="$REF_LEN" '{ sum+=$3; if($3>0) cov++; } \
           END { if (len>0) printf("%.2f\t%d\n", sum/len, cov); else print "0.00\t0"; }')
  fi
  local frac=$(awk -v c="$cov" -v l="$REF_LEN" 'BEGIN{ if(l>0) printf("%.3f", c/l); else print "0.000"; }')

  echo "[$(date +%T)]   Coverage results: Mean Depth = ${mean}x | Covered Bases = ${cov}/${REF_LEN} | Fraction = ${frac}"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$SAMPLE" "$ASSIGNED_TYPE" "$mapped_A" "$mapped_B" "$REF_LEN" "$mean" "$cov" "$frac" >> "$ASSIGNMENT_COVFILE"


  # --- 4) Variant Calling (on assigned BAM) ---
  echo "[$(date +%T)] → Calling variants against $ASSIGNED_REF_NAME"
  bcftools mpileup --fasta-ref "$ASSIGNED_REF" --annotate FORMAT/DP --min-MQ 20 --min-BQ 20 "$ASSIGNED_BAM" \
    | bcftools call --ploidy 1 --variants-only --multiallelic-caller -Oz -o "${OUTDIR}/vcf/${SAMPLE}.assigned.raw.vcf.gz"
  bcftools index "${OUTDIR}/vcf/${SAMPLE}.assigned.raw.vcf.gz"

  # --- 5) Filter Variants ---
  echo "[$(date +%T)] → Filtering variants"
  bcftools filter --exclude 'QUAL<20 || FORMAT/DP<10' \
    "${OUTDIR}/vcf/${SAMPLE}.assigned.raw.vcf.gz" \
    -Oz -o "${OUTDIR}/vcf/${SAMPLE}.assigned.filtered.vcf.gz"
  bcftools index "${OUTDIR}/vcf/${SAMPLE}.assigned.filtered.vcf.gz"

  # --- 6) Generate F-CDS Consensus ---
  echo "[$(date +%T)] → Generating consensus sequence"
  local CONSENSUS_FASTA="${OUTDIR}/consensus/${SAMPLE}.consensus.fasta"
  local CHAIN_FILE="${OUTDIR}/consensus/${SAMPLE}.chain.txt"
  local FILTERED_VCF="${OUTDIR}/vcf/${SAMPLE}.assigned.filtered.vcf.gz"

  # Generate consensus sequence using the assigned reference and filtered variants
  # Removed the problematic '-s "$SAMPLE"' option
  bcftools consensus -f "$ASSIGNED_REF" -H A \
    "$FILTERED_VCF" \
    -c "$CHAIN_FILE" \
    -o "$CONSENSUS_FASTA" \
    || { echo "[$(date +%T)] ERROR: bcftools consensus failed for $SAMPLE"; exit 1; } # Exit if consensus fails

  # Check if consensus file was created and is not empty before renaming header
  if [[ -f "$CONSENSUS_FASTA" && -s "$CONSENSUS_FASTA" ]]; then
    echo "[$(date +%T)]   Consensus file generated successfully. Renaming header."
    local SED_TMP="${CONSENSUS_FASTA}.tmp"
    sed "s/>.*/>${SAMPLE}_F_${ASSIGNED_TYPE}_consensus/" "$CONSENSUS_FASTA" > "$SED_TMP" && mv "$SED_TMP" "$CONSENSUS_FASTA"
  else
    echo "[$(date +%T)] WARNING: Consensus file '$CONSENSUS_FASTA' not created or is empty. Skipping rename and subsequent steps relying on it."
  fi


  # --- 7) Translate Consensus CDS ---
  echo "[$(date +%T)] → Translating consensus CDS"
  local PROTEIN_FASTA="${OUTDIR}/protein/${SAMPLE}.F.protein.fa"
  if [[ -f "$CONSENSUS_FASTA" && -s "$CONSENSUS_FASTA" ]]; then
      transeq -sequence "$CONSENSUS_FASTA" \
              -outseq "$PROTEIN_FASTA" \
              -frame 1 -clean \
              || { echo "[$(date +%T)] ERROR: transeq failed for $SAMPLE. Check '$CONSENSUS_FASTA'."; exit 1; }
      echo "[$(date +%T)]   Protein sequence saved to $PROTEIN_FASTA"
  else
      echo "[$(date +%T)]   Skipping translation because consensus FASTA is missing or empty."
      touch "$PROTEIN_FASTA" # Create empty file as placeholder
  fi

  # --- 8) Annotate Variants using SnpEff ---
  local ANNOTATED_VCF="${OUTDIR}/ann/${SAMPLE}.assigned.ann.vcf"
  local MUTATIONS_TXT="${OUTDIR}/ann/${SAMPLE}.assigned.mutations.txt"
  if [[ $(bcftools view -H "$FILTERED_VCF" | wc -l) -gt 0 ]]; then
      echo "[$(date +%T)] → Annotating variants with SnpEff (DB: $ASSIGNED_SNPEFF_DB)"
      snpEff -noStats -v "$ASSIGNED_SNPEFF_DB" "$FILTERED_VCF" \
        > "$ANNOTATED_VCF" \
        || { echo "[$(date +%T)] WARNING: snpEff annotation failed for $SAMPLE. Check snpEff setup and DB '$ASSIGNED_SNPEFF_DB'."; \
             touch "$ANNOTATED_VCF"; }

      # --- 9) Extract AA changes ---
      if [[ -f "$ANNOTATED_VCF" && -s "$ANNOTATED_VCF" ]]; then
          echo "[$(date +%T)] → Extracting mutations"
          bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%ANN]\n' \
            "$ANNOTATED_VCF" \
            > "$MUTATIONS_TXT" \
            || { echo "[$(date +%T)] WARNING: bcftools query failed for $SAMPLE."; \
                 touch "$MUTATIONS_TXT"; }
          echo "[$(date +%T)]   Annotations saved to $ANNOTATED_VCF"
          echo "[$(date +%T)]   Mutations saved to $MUTATIONS_TXT"
      else
           echo "[$(date +%T)]   Skipping mutation extraction because annotated VCF is missing or empty."
           touch "$MUTATIONS_TXT"
      fi
  else
      echo "[$(date +%T)] → No variants passed filters. Skipping annotation and mutation extraction."
      touch "$ANNOTATED_VCF"
      touch "$MUTATIONS_TXT"
  fi


  echo "[$(date +%T)] === Finished Processing: $SAMPLE (Assigned: $ASSIGNED_TYPE) ==="
  echo "========================================"
}

# Export the function and necessary variables for GNU Parallel
export -f process_sample
export REF_A REF_B REF_A_NAME REF_B_NAME SNPEFF_DB_A SNPEFF_DB_B THREADS OUTDIR ASSIGNMENT_COVFILE

# Determine the number of parallel jobs
TOTAL_CORES=$(nproc)
NUM_JOBS=$(( TOTAL_CORES / 2 ))
if [[ $NUM_JOBS -lt 1 ]]; then NUM_JOBS=1; fi

echo "Starting parallel processing with $NUM_JOBS jobs (each using up to $THREADS threads internally)..."
echo "Input FASTQ files from: filtered_reads/"
echo "Output directory: $OUTDIR"
echo "Log files in: ${OUTDIR}/logs/"

# Find all .fq files in filtered_reads, sort them, and process in parallel
find filtered_reads -name '*.fq' | sort \
  | parallel --eta -j $NUM_JOBS process_sample {} "$REF_A" "$REF_B" "$REF_A_NAME" "$REF_B_NAME" "$SNPEFF_DB_A" "$SNPEFF_DB_B" "$THREADS" "$OUTDIR" "$ASSIGNMENT_COVFILE"

echo "[$(date +%T)] All processing complete."
echo "Assignment and coverage metrics saved to $ASSIGNMENT_COVFILE"

### 10) build a binary mutation matrix of all HGVS.p changes,
###     excluding any sample whose coverage fraction < 0.95 and type is not Undetermined
MATRIX="${OUTDIR}/F_mutation_matrix.tsv"
echo "[$(date +%T)] Building mutation matrix..."

# Use a Python script embedded via HEREDOC for matrix generation
python3 << EOF
import re, glob, pandas as pd
from pathlib import Path
import sys # For stderr

# Define amino acid code mapping (3-letter to 1-letter)
aa3to1 = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C',
    'Glu':'E','Gln':'Q','Gly':'G','His':'H','Ile':'I',
    'Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P',
    'Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V',
    'Ter':'*' # Use '*' for stop codon (termination)
}

# --- Configuration ---
OUTDIR = "${OUTDIR}"
ASSIGNMENT_COVFILE = f"{OUTDIR}/assignment_coverage.tsv"
MATRIX_FILE  = f"{OUTDIR}/F_mutation_matrix.tsv"
MIN_COVERAGE_FRACTION = 0.95 # Minimum coverage fraction threshold

# --- 1) Select samples meeting criteria ---
print(f"Reading assignment/coverage file: {ASSIGNMENT_COVFILE}", file=sys.stderr)
try:
    cov_df = pd.read_csv(ASSIGNMENT_COVFILE, sep='\t')
    cov_df['CoverageFraction'] = pd.to_numeric(cov_df['CoverageFraction'], errors='coerce')
    good_samples_df = cov_df.loc[
        (cov_df['CoverageFraction'] >= MIN_COVERAGE_FRACTION) &
        (cov_df['AssignedType'].isin(['A', 'B']))
    ]
    good_samples = good_samples_df.set_index('Sample')['AssignedType'].to_dict()

except FileNotFoundError:
    print(f"ERROR: Coverage file not found: {ASSIGNMENT_COVFILE}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"ERROR: Could not read or process coverage file: {ASSIGNMENT_COVFILE}", file=sys.stderr)
    print(e, file=sys.stderr)
    sys.exit(1)

if not good_samples:
    print(f"No samples met the coverage threshold ({MIN_COVERAGE_FRACTION}) and had a valid assignment (A/B).", file=sys.stderr)
    print(f"Mutation matrix cannot be generated. Creating empty matrix file: {MATRIX_FILE}", file=sys.stderr)
    with open(MATRIX_FILE, 'w') as f:
         f.write("Mutation\n")
    sys.exit(0)

print(f"Found {len(good_samples)} samples meeting criteria for matrix generation.", file=sys.stderr)

# --- 2) Parse annotations for each good sample ---
pattern = re.compile(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)$')
sample_muts = {}
processed_files_count = 0

print(f"Parsing mutation files from: {OUTDIR}/ann/*.assigned.mutations.txt", file=sys.stderr)
for fn in glob.glob(f"{OUTDIR}/ann/*.assigned.mutations.txt"):
    sample = Path(fn).stem.split('.assigned.mutations')[0]
    if sample not in good_samples:
        continue

    processed_files_count += 1
    assigned_type = good_samples[sample]
    muts_found_in_sample = set()

    try:
        with open(fn) as fh:
            for line_num, line in enumerate(fh, 1):
                if line.startswith('#') or not line.strip():
                    continue
                cols = line.strip().split('\t')
                if len(cols) < 5 or cols[4] in ('.','', None):
                    continue

                annotations = cols[4].split(',')
                for ann_entry in annotations:
                    fields = ann_entry.split('|')
                    if len(fields) > 10:
                        effect = fields[1]
                        hgvs_p = fields[10]

                        if effect in ('missense_variant', 'stop_gained', 'stop_lost', 'start_lost', 'frameshift_variant') and hgvs_p.startswith('p.'):
                            m = pattern.match(hgvs_p)
                            if not m:
                                print(f"Warning: Could not parse HGVS.p format '{hgvs_p}' in sample {sample}, file {fn}, line {line_num}", file=sys.stderr)
                                continue

                            ref_aa3, pos_str, alt_aa3 = m.groups()
                            ref_aa1 = aa3to1.get(ref_aa3)
                            alt_aa1 = aa3to1.get(alt_aa3)

                            if not ref_aa1 or not alt_aa1:
                                print(f"Warning: Unknown amino acid code in '{hgvs_p}' in sample {sample}, file {fn}, line {line_num}", file=sys.stderr)
                                continue
                            if ref_aa1 == alt_aa1:
                                continue

                            pos = int(pos_str)
                            mutation_label = f"{assigned_type}:{ref_aa1}{pos}{alt_aa1}"
                            muts_found_in_sample.add(mutation_label)

    except FileNotFoundError:
        print(f"Warning: Mutation file not found for sample {sample} ({fn}), skipping.", file=sys.stderr)
        continue
    except Exception as e:
        print(f"Warning: Error processing mutation file {fn} for sample {sample}. Error: {e}", file=sys.stderr)
        continue

    if muts_found_in_sample:
        sample_muts[sample] = muts_found_in_sample

print(f"Processed mutation files for {processed_files_count} samples.", file=sys.stderr)

if not sample_muts:
     print(f"No non-synonymous mutations found meeting criteria in any processed sample.", file=sys.stderr)
     print(f"Mutation matrix will be empty. Creating file: {MATRIX_FILE}", file=sys.stderr)
     with open(MATRIX_FILE, 'w') as f:
         f.write("Mutation\t" + "\t".join(sorted(good_samples.keys())) + "\n")
     sys.exit(0)


# --- 3) Build and write the presence/absence matrix ---
print(f"Building mutation matrix...", file=sys.stderr)
all_unique_muts = set()
for muts in sample_muts.values():
    all_unique_muts.update(muts)

def sort_key_mutation(mut_str):
    try:
        type_part, aa_part = mut_str.split(':', 1)
        match = re.search(r'(\d+)', aa_part)
        pos = int(match.group(1)) if match else 0
        return (type_part, pos, aa_part)
    except ValueError:
         print(f"Warning: Could not parse mutation string '{mut_str}' for sorting.", file=sys.stderr)
         return ('', 0, mut_str)

sorted_muts = sorted(list(all_unique_muts), key=sort_key_mutation)
sorted_samples = sorted(sample_muts.keys())

df = pd.DataFrame(0, index=sorted_muts, columns=sorted_samples)
for sample, mutations_in_sample in sample_muts.items():
    for mutation in mutations_in_sample:
        if mutation in df.index:
            df.at[mutation, sample] = 1

df.to_csv(MATRIX_FILE, sep='\t', index_label='Mutation')

print(f"Mutation matrix generated with {len(sorted_muts)} unique mutations and {len(sorted_samples)} samples.", file=sys.stderr)
print(f"Matrix saved to: {MATRIX_FILE}", file=sys.stderr)

# Make sure this EOF is on a line by itself with no leading/trailing whitespace
EOF

# echo complete msg
echo "[$(date +%T)] → Mutation matrix generation complete."
echo "[$(date +%T)] Workflow finished."
