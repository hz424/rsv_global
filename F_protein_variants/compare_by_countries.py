#!/usr/bin/env python3

import pandas as pd
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from datetime import datetime
import sys
import argparse
from collections import defaultdict
import multiprocessing
import subprocess # For calling MAFFT
import tempfile # For temporary files for MAFFT
import io # For reading MAFFT output string as file
import traceback # For detailed error printing

# --- Configuration --- 

# Input directories and files from previous workflow
UAE_WORKFLOW_DIR = Path("workflow_run_AB")
UAE_MATRIX_FILE = UAE_WORKFLOW_DIR / "F_mutation_matrix.tsv"
UAE_ASSIGNMENT_FILE = UAE_WORKFLOW_DIR / "assignment_coverage.tsv"

# Reference sequences used in the previous workflow
REF_A_FASTA = Path("reference/RSV_A_F_cds.fasta")
REF_B_FASTA = Path("reference/RSV_B_F_cds.fasta")

# Public data directory and new files
PUB_DATA_DIR = Path("rsv_F_pub")
PUB_A_METADATA = PUB_DATA_DIR / "metadata_rsvA.tsv"
PUB_A_FASTA = PUB_DATA_DIR / "sequences_rsvA.fasta"
PUB_B_METADATA = PUB_DATA_DIR / "metadata_rsvB.tsv"
PUB_B_FASTA = PUB_DATA_DIR / "sequences_rsvB.fasta"

# Output files
OUTPUT_FREQUENCY_TABLE = Path("rsv_F_mutation_frequency_by_country_mafft_extracted_final.tsv") 
OUTPUT_PARSED_F_PROTEIN_A_FASTA = PUB_DATA_DIR / "parsed_F_proteins_A_mafft_final.fasta" 
OUTPUT_PARSED_F_PROTEIN_B_FASTA = PUB_DATA_DIR / "parsed_F_proteins_B_mafft_final.fasta" 


# Parameters
TARGET_YEARS = [2021,2022,2023,2024]
MIN_F_COVERAGE_PUBLIC = 0.95 
UAE_MIN_DEPTH = 20.0
UAE_MIN_COVERAGE = 0.95

# Length Ratio Thresholds
MIN_EXPECTED_F_GENE_CDS_LENGTH_RATIO = 0.95
MIN_EXPECTED_F_PROTEIN_LENGTH_RATIO = 0.95

NUM_PROCESSES = 50 
PROGRESS_UPDATE_INTERVAL = 100 
MAFFT_TIMEOUT_SECONDS = 120 # Timeout for a single MAFFT alignment

# --- Helper Functions ---

def parse_fasta(fasta_file):
    sequences = {}
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = str(record.seq).upper()
        # print(f"Successfully parsed {len(sequences)} sequences from {fasta_file}", file=sys.stderr) # Can be verbose
    except FileNotFoundError:
        print(f"ERROR: FASTA file not found: {fasta_file}", file=sys.stderr)
        sys.exit(1)
    return sequences

def translate_sequence(nucleotide_seq):
    if not isinstance(nucleotide_seq, str): return ""
    clean_seq = nucleotide_seq.replace("-", "").replace(" ", "")
    if not clean_seq: return ""
    remainder = len(clean_seq) % 3
    if remainder != 0:
        clean_seq = clean_seq[:-remainder]
        if not clean_seq: return ""
    try:
        return str(Seq(clean_seq).translate(to_stop=True))
    except Exception as e:
        return ""

def load_uae_metadata(file_path, columns):
    try:
        df = pd.read_csv(file_path, sep='\t', names=columns, header=0, on_bad_lines='warn')
        print(f"Successfully loaded {len(df)} records from {file_path}", file=sys.stderr)
        return df
    except FileNotFoundError:
        print(f"ERROR: Metadata file not found: {file_path}", file=sys.stderr)
        sys.exit(1)
    return pd.DataFrame()

def load_public_metadata(file_path, required_cols):
    try:
        dtype_spec = {'F_coverage': str}
        df = pd.read_csv(file_path, sep='\t', usecols=required_cols, on_bad_lines='warn', low_memory=False, dtype=dtype_spec)
        print(f"Successfully loaded {len(df)} records from {file_path}", file=sys.stderr)
        return df
    except FileNotFoundError:
        print(f"ERROR: Metadata file not found: {file_path}", file=sys.stderr)
        sys.exit(1)
    except ValueError as ve:
        print(f"ERROR: Problem loading columns from {file_path}. Columns {required_cols}. Error: {ve}", file=sys.stderr)
        sys.exit(1)
    return pd.DataFrame()

def parse_date(date_str):
    if pd.isna(date_str): return pd.NaT
    original_input_date = str(date_str).strip()
    cleaned_date_str = original_input_date.replace('-XX', '-01').replace('/XX', '/01')
    if cleaned_date_str:
        year_part_candidate = re.split(r'[-/]', cleaned_date_str)[0]
        if 'X' in year_part_candidate.upper():
            return pd.NaT
    formats_to_try = ['%Y-%m-%d', '%Y-%m', '%Y', '%m/%d/%y', '%m/%d/%Y', '%Y/%m/%d', '%b-%Y', '%d-%b-%Y']
    for fmt in formats_to_try:
        try:
            return datetime.strptime(cleaned_date_str.split('T')[0], fmt)
        except ValueError:
            continue
    return pd.NaT

def extract_country_simple(location_str):
    if pd.isna(location_str): return "Unknown"
    location_str = str(location_str).strip()
    loc_upper = location_str.upper()
    if "UNITED STATES" in loc_upper or "U.S.A" in loc_upper: return "USA"
    if "UNITED KINGDOM" in loc_upper: return "UK"
    if "SAUDI ARABIA" in loc_upper: return "Saudi Arabia"
    if "SOUTH KOREA" in loc_upper: return "South Korea"
    if "HONG KONG" in loc_upper: return "Hong Kong"
    if "NEW ZEALAND" in loc_upper: return "New Zealand"
    if "SOUTH AFRICA" in loc_upper: return "South Africa"
    if "CENTRAL AFRICAN REPUBLIC" in loc_upper: return "Central African Republic"
    if "EL SALVADOR" in loc_upper: return "El Salvador"
    if "SRI LANKA" in loc_upper: return "Sri Lanka"
    if "BURKINA FASO" in loc_upper: return "Burkina Faso"
    country = re.split(r'\s*[/:\s-]\s*', location_str)[0].strip()
    country_map_abbr = { "USA": "USA", "UK": "UK"}
    if country.upper() in country_map_abbr: return country_map_abbr[country.upper()]
    return country.title() if country else "Unknown"

def find_protein_mutations(sample_protein_seq, ref_protein_seq):
    mutations = set()
    if not isinstance(sample_protein_seq, str) or not isinstance(ref_protein_seq, str): return mutations
    min_len = min(len(sample_protein_seq), len(ref_protein_seq))
    for i in range(min_len):
        ref_aa, sample_aa = ref_protein_seq[i], sample_protein_seq[i]
        if ref_aa in {'-', 'X', 'B', 'Z', 'J', '*'} or sample_aa in {'-', 'X', 'B', 'Z', 'J', '*'}: continue
        if ref_aa != sample_aa: mutations.add(f"{ref_aa}{i+1}{sample_aa}")
    return mutations

def get_pos_from_mutation(mut_str):
    match = re.search(r'\d+', str(mut_str))
    return int(match.group()) if match else 0

def append_protein_to_fasta(filepath, accession_id, protein_sequence_str):
    if not protein_sequence_str: return 
    try:
        record = SeqRecord(Seq(protein_sequence_str), id=accession_id, description="")
        with open(filepath, "a") as f_handle: 
            SeqIO.write([record], f_handle, "fasta")
    except Exception as e:
        print(f"ERROR: Could not append to FASTA file {filepath} for {accession_id}: {e}", file=sys.stderr)

def read_existing_parsed_proteins(fasta_filepath):
    parsed_proteins = {}
    parsed_ids = set()
    if fasta_filepath.exists():
        try:
            for record in SeqIO.parse(fasta_filepath, "fasta"):
                parsed_proteins[record.id] = str(record.seq)
                parsed_ids.add(record.id)
            print(f"Read {len(parsed_ids)} existing parsed proteins from {fasta_filepath}", file=sys.stderr)
        except Exception as e:
            print(f"Warning: Could not properly read existing parsed protein FASTA {fasta_filepath}: {e}.", file=sys.stderr)
            return {}, set() 
    return parsed_proteins, parsed_ids

# --- MAFFT based F-gene extraction (from validated test script) ---
def align_and_extract_with_mafft(ref_cds_seq_record, public_seq_str, public_seq_id_original, is_rc_input=False):
    """
    Aligns ref F-gene to public sequence using MAFFT and extracts F-gene.
    """
    # print(f"\n--- Attempting MAFFT alignment for {public_seq_id_original} (strand: {'RC' if is_rc_input else 'FWD'}) ---") # Verbose
    extracted_f_gene_cds = None
    
    if not public_seq_str:
        # print(f"DEBUG: Public sequence string for {public_seq_id_original} is empty. Skipping MAFFT.", file=sys.stderr)
        return None

    public_seq_id_for_mafft = f"{public_seq_id_original}_pub_target"

    with tempfile.NamedTemporaryFile(mode="w+", delete=True, suffix=".fasta") as tmp_in:
        SeqIO.write([ref_cds_seq_record, SeqRecord(Seq(public_seq_str), id=public_seq_id_for_mafft)], tmp_in.name, "fasta")
        tmp_in.flush()
        # print(f"DEBUG: MAFFT input file '{tmp_in.name}' created with ref '{ref_cds_seq_record.id}' and target '{public_seq_id_for_mafft}'.")

        mafft_cmd = ["mafft", "--localpair", "--quiet", tmp_in.name]
        # print(f"DEBUG: Running MAFFT command: {' '.join(mafft_cmd)}")
        try:
            result = subprocess.run(mafft_cmd, capture_output=True, text=True, check=False, timeout=MAFFT_TIMEOUT_SECONDS, stdin=subprocess.DEVNULL)
            
            # print(f"DEBUG: MAFFT return code for {public_seq_id_original}: {result.returncode}")
            # if result.stdout: print(f"DEBUG: MAFFT stdout (first 100 chars for {public_seq_id_original}):\n{result.stdout[:100]}\n...")
            # if result.stderr: print(f"DEBUG: MAFFT stderr for {public_seq_id_original}:\n{result.stderr}")

            if result.returncode == 0 and result.stdout:
                aligned_seqs_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(io.StringIO(result.stdout), "fasta")}
                
                if ref_cds_seq_record.id in aligned_seqs_dict and public_seq_id_for_mafft in aligned_seqs_dict:
                    aligned_ref_str = aligned_seqs_dict[ref_cds_seq_record.id]
                    aligned_public_str = aligned_seqs_dict[public_seq_id_for_mafft]

                    if len(aligned_ref_str) != len(aligned_public_str):
                        # print(f"DEBUG: Aligned ref and public for {public_seq_id_original} have different lengths from MAFFT.", file=sys.stderr)
                        return None

                    align_len = len(aligned_ref_str)
                    ref_align_start, ref_align_end = -1, -1
                    for i in range(align_len):
                        if aligned_ref_str[i] != '-':
                            if ref_align_start == -1: ref_align_start = i
                            ref_align_end = i
                    
                    if ref_align_start != -1:
                        # print(f"DEBUG: Ref align segment for {public_seq_id_original}: Start {ref_align_start}, End {ref_align_end}")
                        public_segment_corresponding_to_ref = aligned_public_str[ref_align_start : ref_align_end + 1]
                        
                        if is_rc_input:
                            ungapped_public_segment = public_segment_corresponding_to_ref.replace("-", "")
                            extracted_f_gene_cds = str(Seq(ungapped_public_segment).reverse_complement())
                        else:
                            extracted_f_gene_cds = public_segment_corresponding_to_ref.replace("-", "")
                        # print(f"DEBUG: Extracted CDS for {public_seq_id_original} (ungapped, 5'-3'): Length {len(extracted_f_gene_cds)}")
                    # else: print(f"DEBUG: Ref in MAFFT output for {public_seq_id_original} was all gaps.", file=sys.stderr)
                # else: print(f"DEBUG: Expected IDs not in MAFFT output for {public_seq_id_original}.", file=sys.stderr)
            # else: print(f"DEBUG: MAFFT failed or no stdout for {public_seq_id_original}.", file=sys.stderr)
        except subprocess.TimeoutExpired:
            # print(f"DEBUG: MAFFT timed out for {public_seq_id_original}", file=sys.stderr)
            return None
        except FileNotFoundError:
            print(f"CRITICAL ERROR: MAFFT command not found during alignment of {public_seq_id_original}. Ensure MAFFT is installed and in PATH.", file=sys.stderr)
            # This error should ideally be caught once at the start, but good to have a fallback.
            return None 
        except Exception as e:
            # print(f"DEBUG: Error with MAFFT for {public_seq_id_original}: {type(e).__name__}: {e}", file=sys.stderr)
            return None
            
    return extracted_f_gene_cds

# --- Worker function for parallel processing ---
def process_public_sequence(args_tuple):
    """
    Worker function to process a single public sequence using MAFFT.
    Extracts F-gene, translates, and finds mutations.
    """
    accession_id, full_nt_public_str, country, rsv_type, \
    ref_cds_seq_record_obj, ref_prot_for_comparison, \
    min_expected_f_gene_cds_len_abs, min_expected_f_prot_len_abs = args_tuple

    prot_seq_str = None
    public_f_gene_cds_extracted = None
    
    # Attempt 1: Forward strand
    extracted_fwd = align_and_extract_with_mafft(ref_cds_seq_record_obj, full_nt_public_str, accession_id, is_rc_input=False)
    
    if extracted_fwd and len(extracted_fwd) >= min_expected_f_gene_cds_len_abs:
        public_f_gene_cds_extracted = extracted_fwd
    else:
        # Attempt 2: Reverse complement strand if forward failed or was too short
        public_seq_obj_for_rc = Seq(full_nt_public_str)
        public_seq_rc_str = str(public_seq_obj_for_rc.reverse_complement())
        extracted_rev = align_and_extract_with_mafft(ref_cds_seq_record_obj, public_seq_rc_str, accession_id, is_rc_input=True)
        
        if extracted_rev and len(extracted_rev) >= min_expected_f_gene_cds_len_abs:
            public_f_gene_cds_extracted = extracted_rev

    if not public_f_gene_cds_extracted:
        return accession_id, rsv_type, country, set(), "no_mafft_alignment_or_short_cds", prot_seq_str

    prot_seq_str = translate_sequence(public_f_gene_cds_extracted)
    if not prot_seq_str:
        return accession_id, rsv_type, country, set(), "translation_failed", None
    
    if len(prot_seq_str) < min_expected_f_prot_len_abs:
        return accession_id, rsv_type, country, set(), "short_f_protein", None

    muts = find_protein_mutations(prot_seq_str, ref_prot_for_comparison)
    return accession_id, rsv_type, country, {f"{rsv_type}:{m}" for m in muts}, "processed_mafft", prot_seq_str


# --- Main Script ---
def main():
    print("Starting RSV F Mutation Frequency Comparison Analysis (MAFFT - Refined Extraction)...")

    # MAFFT Availability Check (from test script)
    mafft_ok = False
    try:
        process = subprocess.run(["mafft", "--version"], capture_output=True, text=True, timeout=10, stdin=subprocess.DEVNULL)
        if process.returncode == 0 and process.stderr:
            # print(f"DEBUG: MAFFT version check successful. STDERR:\n{process.stderr.strip()}", file=sys.stderr)
            if "mafft" in process.stderr.lower() or re.search(r"v\d+\.\d+", process.stderr):
                 mafft_ok = True
            else: # RC=0, stderr exists but doesn't look like version info
                 mafft_ok = True # Assume ok if it ran and produced stderr
        elif process.returncode == 0 and process.stdout: # Fallback: check stdout
            # print(f"DEBUG: MAFFT version check successful (RC=0, output to STDOUT). STDOUT:\n{process.stdout.strip()}", file=sys.stderr)
            if "mafft" in process.stdout.lower() or re.search(r"v\d+\.\d+", process.stdout):
                 mafft_ok = True
            else: # RC=0, stdout exists but no clear version
                 mafft_ok = True # Assume ok
        else: # Non-zero RC or no output
            print("CRITICAL ERROR: `mafft --version` failed or produced no recognizable output.", file=sys.stderr)
            print(f"  Return Code: {process.returncode}\n  STDERR: {process.stderr.strip()}\n  STDOUT: {process.stdout.strip()}", file=sys.stderr)
            sys.exit(1)
    except FileNotFoundError:
        print("CRITICAL ERROR: MAFFT command not found. Please ensure MAFFT is installed and in your system PATH.", file=sys.stderr)
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print("CRITICAL ERROR: MAFFT check (`mafft --version`) timed out.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"CRITICAL ERROR: An unexpected error occurred during MAFFT check: {type(e).__name__}: {e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)
    
    if not mafft_ok: # Should have exited above if not OK, but as a safeguard.
        print("CRITICAL ERROR: MAFFT not found or check failed. Exiting.", file=sys.stderr)
        sys.exit(1)
    print("MAFFT availability check successful.", file=sys.stderr)


    # 1. Load UAE Data and Apply Filters
    print("\n--- Loading and Filtering UAE Data ---")
    try:
        uae_matrix_full_df = pd.read_csv(UAE_MATRIX_FILE, sep='\t', index_col='Mutation')
        print(f"Loaded full UAE mutation matrix: {uae_matrix_full_df.shape[0]} mutations, {uae_matrix_full_df.shape[1]} samples", file=sys.stderr)
    except FileNotFoundError:
        print(f"ERROR: UAE Mutation Matrix file not found: {UAE_MATRIX_FILE}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Could not read UAE mutation matrix file {UAE_MATRIX_FILE}: {e}", file=sys.stderr)
        sys.exit(1)
        
    uae_assignment_df = load_uae_metadata(UAE_ASSIGNMENT_FILE, ['Sample', 'AssignedType', 'MappedReads_A', 'MappedReads_B', 'ReferenceLength', 'MeanDepth', 'CoveredBases', 'CoverageFraction'])
    uae_assignment_df['MeanDepth'] = pd.to_numeric(uae_assignment_df['MeanDepth'], errors='coerce')
    uae_assignment_df['CoverageFraction'] = pd.to_numeric(uae_assignment_df['CoverageFraction'], errors='coerce')

    print(f"Filtering UAE samples: MeanDepth > {UAE_MIN_DEPTH} AND CoverageFraction > {UAE_MIN_COVERAGE}", file=sys.stderr)
    initial_uae_samples = len(uae_assignment_df)
    uae_assignment_filtered_df = uae_assignment_df[
        (uae_assignment_df['MeanDepth'] > UAE_MIN_DEPTH) &
        (uae_assignment_df['CoverageFraction'] > UAE_MIN_COVERAGE)
    ].copy()
    print(f"Retained {len(uae_assignment_filtered_df)} UAE samples out of {initial_uae_samples} after filtering.", file=sys.stderr)

    filtered_uae_sample_names = uae_assignment_filtered_df['Sample'].tolist()
    valid_filtered_uae_samples = [s for s in filtered_uae_sample_names if s in uae_matrix_full_df.columns]
    if not valid_filtered_uae_samples and filtered_uae_sample_names: 
        uae_matrix_df = pd.DataFrame() 
    elif not filtered_uae_sample_names:
        uae_matrix_df = pd.DataFrame()
    else:
        uae_matrix_df = uae_matrix_full_df[valid_filtered_uae_samples]
    print(f"Filtered UAE mutation matrix to {uae_matrix_df.shape[1]} samples.", file=sys.stderr)
    uae_assignment_map = uae_assignment_filtered_df.set_index('Sample')['AssignedType'].to_dict()


    # 2. Load Reference F-gene Protein and Nucleotide Sequences
    print("\n--- Loading Reference F-gene Sequences and Translating ---")
    ref_a_nt_map = parse_fasta(REF_A_FASTA) 
    ref_b_nt_map = parse_fasta(REF_B_FASTA)
    if not ref_a_nt_map or not ref_b_nt_map: sys.exit(1)

    ref_a_id, ref_b_id = list(ref_a_nt_map.keys())[0], list(ref_b_nt_map.keys())[0]
    ref_a_nt_cds_str, ref_b_nt_cds_str = ref_a_nt_map[ref_a_id], ref_b_nt_map[ref_b_id]
    
    ref_a_cds_seq_record = SeqRecord(Seq(ref_a_nt_cds_str), id=ref_a_id, name=ref_a_id, description="Ref_A_F_CDS")
    ref_b_cds_seq_record = SeqRecord(Seq(ref_b_nt_cds_str), id=ref_b_id, name=ref_b_id, description="Ref_B_F_CDS")

    ref_a_prot, ref_b_prot = translate_sequence(ref_a_nt_cds_str), translate_sequence(ref_b_nt_cds_str)
    if not ref_a_prot or not ref_b_prot:
        print("ERROR: Failed to translate one or both reference F-gene sequences. Exiting.", file=sys.stderr)
        sys.exit(1)
    print(f"Reference A F-protein Length: {len(ref_a_prot)}, CDS Length: {len(ref_a_nt_cds_str)}", file=sys.stderr)
    print(f"Reference B F-protein Length: {len(ref_b_prot)}, CDS Length: {len(ref_b_nt_cds_str)}", file=sys.stderr)
    
    min_expected_f_gene_cds_len_a_abs = len(ref_a_nt_cds_str) * MIN_EXPECTED_F_GENE_CDS_LENGTH_RATIO
    min_expected_f_gene_cds_len_b_abs = len(ref_b_nt_cds_str) * MIN_EXPECTED_F_GENE_CDS_LENGTH_RATIO
    min_expected_f_prot_len_a_abs = len(ref_a_prot) * MIN_EXPECTED_F_PROTEIN_LENGTH_RATIO
    min_expected_f_prot_len_b_abs = len(ref_b_prot) * MIN_EXPECTED_F_PROTEIN_LENGTH_RATIO


    # 3. Load and Filter Public Data Metadata
    print("\n--- Loading and Filtering Public Data Metadata ---")
    public_metadata_cols = ['accession', 'date', 'country', 'F_coverage']
    
    pub_a_meta_df = load_public_metadata(PUB_A_METADATA, public_metadata_cols)
    pub_a_meta_recent_df = pd.DataFrame() 
    if not pub_a_meta_df.empty:
        pub_a_meta_df['F_coverage'] = pd.to_numeric(pub_a_meta_df['F_coverage'], errors='coerce').fillna(0.0)
        pub_a_meta_df = pub_a_meta_df[pub_a_meta_df['F_coverage'] >= MIN_F_COVERAGE_PUBLIC].copy()
        pub_a_meta_df['parsed_date'] = pub_a_meta_df['date'].apply(parse_date)
        pub_a_meta_df = pub_a_meta_df.dropna(subset=['parsed_date'])
        pub_a_meta_recent_df = pub_a_meta_df[pub_a_meta_df['parsed_date'].dt.year.isin(TARGET_YEARS)].copy()
        if not pub_a_meta_recent_df.empty:
            pub_a_meta_recent_df['country_clean'] = pub_a_meta_recent_df['country'].apply(extract_country_simple)
        else:
            pub_a_meta_recent_df['country_clean'] = pd.Series(dtype='str')
    print(f"RSV-A Metadata: Found {len(pub_a_meta_recent_df)} entries for years {TARGET_YEARS} after filters.", file=sys.stderr)

    pub_b_meta_df = load_public_metadata(PUB_B_METADATA, public_metadata_cols)
    pub_b_meta_recent_df = pd.DataFrame() 
    if not pub_b_meta_df.empty:
        pub_b_meta_df['F_coverage'] = pd.to_numeric(pub_b_meta_df['F_coverage'], errors='coerce').fillna(0.0)
        pub_b_meta_df = pub_b_meta_df[pub_b_meta_df['F_coverage'] >= MIN_F_COVERAGE_PUBLIC].copy()
        pub_b_meta_df['parsed_date'] = pub_b_meta_df['date'].apply(parse_date)
        pub_b_meta_df = pub_b_meta_df.dropna(subset=['parsed_date'])
        pub_b_meta_recent_df = pub_b_meta_df[pub_b_meta_df['parsed_date'].dt.year.isin(TARGET_YEARS)].copy()
        if not pub_b_meta_recent_df.empty:
            pub_b_meta_recent_df['country_clean'] = pub_b_meta_recent_df['country'].apply(extract_country_simple)
        else:
            pub_b_meta_recent_df['country_clean'] = pd.Series(dtype='str')
    print(f"RSV-B Metadata: Found {len(pub_b_meta_recent_df)} entries for years {TARGET_YEARS} after filters.", file=sys.stderr)


    # 4. Load Public Sequences, Extract F-gene, and Identify Mutations (Parallel with Resume)
    print("\n--- Processing Public Sequences (Parallel with Resume): Extracting F-gene and Identifying Mutations ---")
    pub_mutations = defaultdict(set)
    sample_to_type = {}
    sample_to_country = {}
    
    existing_parsed_a_proteins, already_parsed_a_ids = read_existing_parsed_proteins(OUTPUT_PARSED_F_PROTEIN_A_FASTA)
    existing_parsed_b_proteins, already_parsed_b_ids = read_existing_parsed_proteins(OUTPUT_PARSED_F_PROTEIN_B_FASTA)
    
    print("Processing mutations for already parsed RSV-A proteins...", file=sys.stderr)
    processed_from_existing_a = 0
    for acc_id, prot_seq in existing_parsed_a_proteins.items():
        if not pub_a_meta_recent_df.empty and acc_id in pub_a_meta_recent_df['accession'].values:
            if len(prot_seq) < min_expected_f_prot_len_a_abs:
                continue 

            country_val = pub_a_meta_recent_df.loc[pub_a_meta_recent_df['accession'] == acc_id, 'country_clean'].iloc[0]
            muts = find_protein_mutations(prot_seq, ref_a_prot)
            pub_mutations[acc_id].update({f"A:{m}" for m in muts})
            sample_to_type[acc_id] = 'A'
            sample_to_country[acc_id] = country_val
            processed_from_existing_a += 1
    print(f"  Finished processing {processed_from_existing_a} relevant already parsed RSV-A proteins for mutations.", file=sys.stderr)

    processed_from_existing_b = 0
    print("Processing mutations for already parsed RSV-B proteins...", file=sys.stderr)
    for acc_id, prot_seq in existing_parsed_b_proteins.items():
        if not pub_b_meta_recent_df.empty and acc_id in pub_b_meta_recent_df['accession'].values:
            if len(prot_seq) < min_expected_f_prot_len_b_abs:
                continue

            country_val = pub_b_meta_recent_df.loc[pub_b_meta_recent_df['accession'] == acc_id, 'country_clean'].iloc[0]
            muts = find_protein_mutations(prot_seq, ref_b_prot)
            pub_mutations[acc_id].update({f"B:{m}" for m in muts})
            sample_to_type[acc_id] = 'B'
            sample_to_country[acc_id] = country_val
            processed_from_existing_b +=1
    print(f"  Finished processing {processed_from_existing_b} relevant already parsed RSV-B proteins for mutations.", file=sys.stderr)


    tasks_a, tasks_b = [], []

    if not pub_a_meta_recent_df.empty:
        print("\nPreparing NEW tasks for Public RSV-A sequences...")
        pub_a_full_seqs_nt = parse_fasta(PUB_A_FASTA)
        pub_a_meta_recent_unique = pub_a_meta_recent_df.drop_duplicates(subset=['accession'], keep='first')
        if 'country_clean' not in pub_a_meta_recent_unique.columns:
             pub_a_meta_recent_unique['country_clean'] = "Unknown" 

        pub_a_country_map = pub_a_meta_recent_unique.set_index('accession')['country_clean'].to_dict()
        
        for accession_id in pub_a_meta_recent_unique['accession']:
            if accession_id not in already_parsed_a_ids and accession_id in pub_a_full_seqs_nt: 
                country = pub_a_country_map.get(accession_id, "Unknown")
                tasks_a.append((accession_id, pub_a_full_seqs_nt[accession_id], country, 'A',
                                ref_a_cds_seq_record, ref_a_prot, min_expected_f_gene_cds_len_a_abs, min_expected_f_prot_len_a_abs))
        print(f"  Prepared {len(tasks_a)} new tasks for RSV-A.", file=sys.stderr)


    if not pub_b_meta_recent_df.empty:
        print("\nPreparing NEW tasks for Public RSV-B sequences...")
        pub_b_full_seqs_nt = parse_fasta(PUB_B_FASTA)
        pub_b_meta_recent_unique = pub_b_meta_recent_df.drop_duplicates(subset=['accession'], keep='first')
        if 'country_clean' not in pub_b_meta_recent_unique.columns:
             pub_b_meta_recent_unique['country_clean'] = "Unknown"

        pub_b_country_map = pub_b_meta_recent_unique.set_index('accession')['country_clean'].to_dict()
        for accession_id in pub_b_meta_recent_unique['accession']:
            if accession_id not in already_parsed_b_ids and accession_id in pub_b_full_seqs_nt: 
                country = pub_b_country_map.get(accession_id, "Unknown")
                tasks_b.append((accession_id, pub_b_full_seqs_nt[accession_id], country, 'B',
                                ref_b_cds_seq_record, ref_b_prot, min_expected_f_gene_cds_len_b_abs, min_expected_f_prot_len_b_abs))
        print(f"  Prepared {len(tasks_b)} new tasks for RSV-B.", file=sys.stderr)


    all_tasks = tasks_a + tasks_b
    total_tasks_to_process = len(all_tasks)
    
    processed_count_in_pool = 0
    skipped_stats = defaultdict(int)

    if all_tasks:
        print(f"\nProcessing {total_tasks_to_process} NEW public sequences using {NUM_PROCESSES} cores...")
        with multiprocessing.Pool(processes=NUM_PROCESSES) as pool:
            for result in pool.imap_unordered(process_public_sequence, all_tasks):
                processed_count_in_pool += 1
                acc_id, rsv_type, country_val, mutations_set, status, prot_seq = result
                
                if status == "processed_mafft": 
                    pub_mutations[acc_id].update(mutations_set)
                    sample_to_type[acc_id] = rsv_type
                    sample_to_country[acc_id] = country_val
                    if rsv_type == 'A' and prot_seq: 
                        append_protein_to_fasta(OUTPUT_PARSED_F_PROTEIN_A_FASTA, acc_id, prot_seq)
                    elif rsv_type == 'B' and prot_seq:
                        append_protein_to_fasta(OUTPUT_PARSED_F_PROTEIN_B_FASTA, acc_id, prot_seq)
                else:
                    skipped_stats[status] += 1
                
                if processed_count_in_pool % PROGRESS_UPDATE_INTERVAL == 0 or processed_count_in_pool == total_tasks_to_process:
                    print(f"  Progress (new tasks): Processed {processed_count_in_pool}/{total_tasks_to_process} sequences...", file=sys.stderr, flush=True)
        
        print("Finished parallel processing of NEW public sequences.")
    else:
        print("No NEW public sequences to process via the pool.")

    total_successfully_processed_public = len(sample_to_type) 

    print(f"\nTotal successfully analyzed F-genes (including previously parsed): {total_successfully_processed_public}", file=sys.stderr)
    for reason, count_val in skipped_stats.items(): 
        print(f"  Skipped {count_val} NEW public sequences (during pool processing) due to: {reason}", file=sys.stderr)
    
    if total_successfully_processed_public > 0:
        unique_countries_assigned = sorted(list(set(c for c in sample_to_country.values() if c != "Unknown")))
        print(f"Unique known countries assigned to all processed public sequences: {unique_countries_assigned}", file=sys.stderr)

    # 5. Calculate UAE Mutation Frequencies
    print("\n--- Calculating UAE Mutation Frequencies (Filtered Samples) ---")
    uae_freq_list = []
    total_uae_a_filtered = sum(1 for st in uae_assignment_map.values() if st == 'A')
    total_uae_b_filtered = sum(1 for st in uae_assignment_map.values() if st == 'B')
    if not uae_matrix_df.empty:
        for mutation, row in uae_matrix_df.iterrows():
            mutation_type = mutation.split(':')[0]
            count = row.sum()
            if mutation_type == 'A' and total_uae_a_filtered > 0:
                uae_freq_list.append({'Mutation': mutation, 'Type': 'A', 'Country': 'UAE', 'Frequency': count / total_uae_a_filtered, 'Count': int(count), 'Total_Samples': total_uae_a_filtered})
            elif mutation_type == 'B' and total_uae_b_filtered > 0:
                uae_freq_list.append({'Mutation': mutation, 'Type': 'B', 'Country': 'UAE', 'Frequency': count / total_uae_b_filtered, 'Count': int(count), 'Total_Samples': total_uae_b_filtered})
    uae_freq_df = pd.DataFrame(uae_freq_list)
    print(f"Calculated {len(uae_freq_df)} UAE mutation frequency entries.", file=sys.stderr)


    # 6. Calculate Public Mutation Frequencies by Country
    print("\n--- Calculating Public Mutation Frequencies by Country ---")
    pub_freq_list = []
    if total_successfully_processed_public > 0: 
        country_type_samples = defaultdict(lambda: defaultdict(list))
        for sample_id, s_type_val in sample_to_type.items():
            country_val = sample_to_country.get(sample_id, "Unknown")
            country_type_samples[country_val][s_type_val].append(sample_id)
        
        for country_val, type_dict in country_type_samples.items():
            for s_type_val, samples_in_group in type_dict.items():
                total_samples_in_group = len(samples_in_group)
                if total_samples_in_group == 0: continue
                mutations_in_group_counts = defaultdict(int)
                for sample_id in samples_in_group:
                    for mut_str in pub_mutations.get(sample_id, set()):
                        if mut_str.startswith(f"{s_type_val}:"):
                            mutations_in_group_counts[mut_str] += 1
                for mut_str, count_val in mutations_in_group_counts.items():
                    pub_freq_list.append({'Mutation': mut_str, 'Type': s_type_val, 'Country': country_val, 'Frequency': count_val / total_samples_in_group, 'Count': count_val, 'Total_Samples': total_samples_in_group})
    pub_freq_df = pd.DataFrame(pub_freq_list)
    print(f"Calculated {len(pub_freq_df)} public mutation frequency entries.", file=sys.stderr)


    # 7. Combine and Save Results
    print("\n--- Combining Frequencies and Saving Results ---")
    combined_df = pd.concat([uae_freq_df, pub_freq_df], ignore_index=True)
    if combined_df.empty:
        print("Combined table is empty. No output file will be generated.", file=sys.stderr)
    else:
        combined_df['Count'] = combined_df['Count'].astype(int)
        combined_df['Total_Samples'] = combined_df['Total_Samples'].astype(int)
        combined_df['Frequency'] = combined_df['Frequency'].astype(float)
        combined_df['Pos'] = combined_df['Mutation'].apply(get_pos_from_mutation)
        final_sorted_df = combined_df.sort_values(by=['Type', 'Pos', 'Country', 'Mutation']).drop(columns=['Pos'])
        if not final_sorted_df.empty:
             print(f"Unique countries in the final output table: {sorted(final_sorted_df['Country'].unique())}", file=sys.stderr)
        else:
            print("Final sorted DataFrame is empty, no frequency table to save.", file=sys.stderr)

        try:
            if not final_sorted_df.empty:
                final_sorted_df.to_csv(OUTPUT_FREQUENCY_TABLE, sep='\t', index=False, float_format='%.4f')
                print(f"\nFrequency table saved successfully to: {OUTPUT_FREQUENCY_TABLE}", file=sys.stderr)
        except Exception as e:
            print(f"\nERROR: Could not save output file {OUTPUT_FREQUENCY_TABLE}: {e}", file=sys.stderr)
            sys.exit(1)
            
    print("\nAnalysis finished.")


if __name__ == '__main__':
    main()
