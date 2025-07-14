# analyze_openmm.py
# FINAL version with corrected RMSD time axis and annotations.
import sys
import os
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt
import numpy as np

def analyze_trajectory(topology_file, trajectory_file, reference_pdb, label):
    # Check if input files exist before trying to load them
    if not all(os.path.exists(f) for f in [topology_file, trajectory_file, reference_pdb]):
        print(f"Error: One or more input files for '{label}' not found. Skipping analysis.")
        return None, None, None

    # Load the trajectory and the stable reference structure
    u = mda.Universe(topology_file, trajectory_file)
    ref = mda.Universe(reference_pdb)

    # --- Selections ---
    protein_complex = u.select_atoms("protein")
    f_protein_chain_A_traj = u.select_atoms("protein and segid A")
    f_protein_chain_A_ref = ref.select_atoms("protein and segid A")
    
    # --- Correct RMSD Calculation ---
    # We align the trajectory to the reference PDB and calculate RMSD against it.
    R_complex = rms.RMSD(protein_complex, ref.select_atoms("protein"), select='backbone', groupselections=['backbone and name CA', 'backbone and name CA'])
    R_complex.run()
    rmsd_data = R_complex.results.rmsd
    
    # --- RMSF Calculation on F-protein from the trajectory ---
    calphas_chain_A = f_protein_chain_A_traj.select_atoms('name CA')
    R_fluct = rms.RMSF(calphas_chain_A).run()
    rmsf_data = R_fluct.results.rmsf
    
    # --- Get PDB Residue IDs ---
    pdb_resids = f_protein_chain_A_ref.select_atoms('name CA').resids
    
    # --- Data preparation for plotting with gaps ---
    full_res_range = np.arange(pdb_resids.min(), pdb_resids.max() + 1)
    gapped_rmsf = np.full(full_res_range.shape, np.nan)
    res_map = {resid: i for i, resid in enumerate(full_res_range)}
    for resid, rsmf_val in zip(pdb_resids, rmsf_data):
        if resid in res_map:
            gapped_rmsf[res_map[resid]] = rsmf_val
            
    return rmsd_data, gapped_rmsf, full_res_range

print("Analyzing Wild-Type...")
wt_rmsd, wt_gapped_rmsf, wt_full_range = analyze_trajectory(
    './wt_minimized.pdb', './wt_trajectory.dcd', './5udd_complex_wt_clean.pdb', 'Wild-Type'
)

print("\nAnalyzing Poly-Mutant...")
mut_rmsd, mut_gapped_rmsf, mut_full_range = analyze_trajectory(
    './mut_minimized.pdb', './mut_trajectory.dcd', './5udd_polymutant.pdb', 'Poly-Mutant'
)

if all(data is not None for data in [wt_rmsd, mut_rmsd, wt_gapped_rmsf, mut_gapped_rmsf]):
    print("\nAll analyses completed successfully. Generating plots...")
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 14), gridspec_kw={'height_ratios': [1, 1]})

    # --- Corrected time axis calculation ---
    wt_time_raw = wt_rmsd[:, 1]
    mut_time_raw = mut_rmsd[:, 1]
    # Subtract the first time value to start the axis at 0, then convert to ns
    wt_time_ns = (wt_time_raw - wt_time_raw[0]) / 1000
    mut_time_ns = (mut_time_raw - mut_time_raw[0]) / 1000
    
    # Plot RMSD
    ax1.plot(wt_time_ns, wt_rmsd[:, 2], label='Wild-Type', color='blue', alpha=0.8, linewidth=1.5) 
    ax1.plot(mut_time_ns, mut_rmsd[:, 2], label='Poly-Mutant', color='red', alpha=0.8, linewidth=1.5)
    ax1.set_title('Backbone RMSD of the Complex', fontsize=24)
    ax1.set_xlabel('Time (ns)', fontsize=20)
    ax1.set_ylabel('RMSD (Å)', fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=18)
    ax1.legend(fontsize=20)
    ax1.set_ylim(0) 

    # Plot RMSF
    ax2.plot(wt_full_range, wt_gapped_rmsf, label='Wild-Type', color='blue', alpha=0.7, linewidth=1.5)
    ax2.plot(mut_full_range, mut_gapped_rmsf, label='Poly-Mutant', color='red', alpha=0.7, linewidth=1.5)
    ax2.set_title('C-Alpha RMSF of RSV F-Protein (Chain A)', fontsize=24)
    ax2.set_xlabel('Residue ID (PDB Numbering)', fontsize=20)
    ax2.set_ylabel('RMSF (Å)', fontsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=18)
    ax2.legend(fontsize=20)
    ax2.set_ylim(0)
    
    # --- Automatically label the mutation peaks ---
    mutation_sites = {'S173P': 173, 'I206M': 206, 'Q209R': 209, 'S211N': 211, 'V220I': 220}
    label_offsets = {'S173P': (0, 1.5), 'I206M': (-25, 3), 'Q209R': (0, 1.5), 'S211N': (25, 3), 'V220I': (0, 1.5)}
    
    for label, res_id in mutation_sites.items():
        try:
            y_coord = mut_gapped_rmsf[np.where(mut_full_range == res_id)[0][0]]
            if np.isnan(y_coord): continue
            x_offset, y_offset = label_offsets.get(label, (0, 1))
            ax2.annotate(label, 
                         xy=(res_id, y_coord), 
                         xytext=(res_id + x_offset, y_coord + y_offset),
                         arrowprops=dict(facecolor='black', shrink=0.05, width=0.5, headwidth=8, connectionstyle="arc3,rad=0.1"),
                         ha='center', va='bottom', fontsize=16, 
                         bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.7))
        except IndexError:
            print(f"Warning: Residue ID {res_id} for mutation {label} not found in trajectory data.")


    plt.tight_layout()
    output_filename = 'md_analysis_comparison_openmm.svg'
    plt.savefig(output_filename, format='svg')
    print(f"Analysis plots saved to {output_filename}")
else:
    print("\nAnalysis could not be completed due to missing simulation files. Exiting.")
