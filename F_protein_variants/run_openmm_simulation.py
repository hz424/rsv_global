# run_openmm_simulation.py
#
# A complete script to run an MD simulation using OpenMM.
# This version uses a robust checkpointing system to resume interrupted runs
# and has been updated for longer, more stable simulations.

import sys
import os
from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer

def run_simulation(pdb_filename, output_prefix):
    """
    Runs a full MD simulation workflow for a given PDB file,
    resuming from a checkpoint if possible.

    This revised function runs a much longer simulation (50 ns)
    to ensure the system reaches a stable equilibrium (plateau).
    """
    print(f"--- Starting Simulation for {pdb_filename} ---")

    # --- Simulation Parameters ---
    # These are now defined at the top for easy modification.
    time_step = 0.002 * picoseconds
    temperature = 300 * kelvin
    equilibration_steps = 250000  # 500 ps equilibration
    production_steps = 25000000  # 50 ns production run (increased from 1 ns)
    reporting_interval = 5000     # Report state data every 10 ps
    checkpoint_interval = 25000  # Save a checkpoint every 50 ps

    equilibration_time = (equilibration_steps * time_step).value_in_unit(nanoseconds)
    production_time = (production_steps * time_step).value_in_unit(nanoseconds)
    print(f"Simulation Configuration: {equilibration_time:.2f} ns Equilibration, {production_time:.0f} ns Production")


    # --- File Existence Checks ---
    final_trajectory_file = f'{output_prefix}_trajectory.dcd'
    checkpoint_file = f'{output_prefix}_checkpoint.chk'
    minimized_pdb_file = f'{output_prefix}_minimized.pdb'

    if os.path.exists(final_trajectory_file):
        print(f"Final output file '{final_trajectory_file}' already exists. Skipping simulation.")
        return

    # --- 1. System Setup ---
    print("\n1. Setting up system...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # If a checkpoint or minimized PDB exists, we don't need to re-run PDBFixer.
    # The topology can be loaded from the minimized PDB.
    if not os.path.exists(checkpoint_file) and not os.path.exists(minimized_pdb_file):
        print("   > No checkpoint or minimized file found. Building system from original PDB...")
        fixer = PDBFixer(filename=pdb_filename)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(True)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)

        modeller = Modeller(fixer.topology, fixer.positions)

        print("   > Adding solvent...")
        # Using a larger padding distance for the solvent box for more stability.
        modeller.addSolvent(forcefield, padding=1.2*nanometers, model='tip3p')

        topology = modeller.topology
        positions = modeller.positions
        # Write the initial PDB with solvent for inspection
        PDBFile.writeFile(topology, positions, open(f'{output_prefix}_initial_solvated.pdb', 'w'))

    else:
        print(f"   > Checkpoint/minimized PDB found. System topology will be loaded from '{minimized_pdb_file}'.")
        # Ensure the minimized file actually exists before trying to load it.
        if not os.path.exists(minimized_pdb_file):
             print(f"ERROR: Checkpoint exists but '{minimized_pdb_file}' is missing. Cannot resume.")
             sys.exit(1)
        pdb = PDBFile(minimized_pdb_file)
        topology = pdb.topology
        positions = pdb.positions

    system = forcefield.createSystem(topology, nonbondedMethod=PME,
                                     nonbondedCutoff=1.0*nanometers, constraints=HBonds)
    system.addForce(MonteCarloBarostat(1*bar, temperature, 25)) # Add barostat for NPT ensemble

    # --- 2. Simulation Platform Setup ---
    print("\n2. Setting up simulation platform...")
    integrator = LangevinMiddleIntegrator(temperature, 1/picosecond, time_step)
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'mixed'}
        print("   > Found and using CUDA platform (NVIDIA GPU) with mixed precision.")
        simulation = Simulation(topology, system, integrator, platform, properties)
    except Exception:
        try:
            platform = Platform.getPlatformByName('OpenCL')
            properties = {'Precision': 'mixed'}
            print("   > Found and using OpenCL platform (AMD/Intel GPU) with mixed precision.")
            simulation = Simulation(topology, system, integrator, platform, properties)
        except Exception:
            platform = Platform.getPlatformByName('CPU')
            print("   > WARNING: No GPU platform found. Using CPU. Simulation will be VERY slow.")
            simulation = Simulation(topology, system, integrator, platform)


    simulation.context.setPositions(positions)

    # --- Steps 3 & 4 are only run if not resuming ---
    if not os.path.exists(checkpoint_file):
        # --- 3. Energy Minimization ---
        print("\n3. Minimizing energy...")
        simulation.minimizeEnergy()
        min_positions = simulation.context.getState(getPositions=True).getPositions()
        # Save the minimized PDB file. This is crucial for resuming.
        PDBFile.writeFile(simulation.topology, min_positions, open(minimized_pdb_file, 'w'))
        print(f"   > Minimized coordinates saved to '{minimized_pdb_file}'.")

        # --- 4. NPT Equilibration ---
        print(f"\n4. Equilibrating system for {equilibration_time:.2f} ns...")
        simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval, step=True,
                                                       potentialEnergy=True, temperature=True, progress=True,
                                                       remainingTime=True, speed=True, totalSteps=equilibration_steps, separator='\t'))
        simulation.step(equilibration_steps)

        # Save the first checkpoint after equilibration
        simulation.saveCheckpoint(checkpoint_file)
        print(f"\n   > Equilibration complete. Checkpoint saved to '{checkpoint_file}'.")

    # --- 5. Production Simulation ---
    print(f"\n5. Starting {production_time:.0f} ns Production MD...")
    # Load state from checkpoint if it exists. This ensures we start where we left off.
    if os.path.exists(checkpoint_file):
        print(f"   > Loading state from checkpoint: '{checkpoint_file}'")
        simulation.loadCheckpoint(checkpoint_file)

    # Clear old reporters and set up new ones for the production run
    simulation.reporters.clear()
    simulation.reporters.append(DCDReporter(final_trajectory_file, reporting_interval))
    simulation.reporters.append(StateDataReporter(f'{output_prefix}_log.txt', reporting_interval, step=True,
                                                   time=True, potentialEnergy=True, temperature=True,
                                                   volume=True, speed=True))
    simulation.reporters.append(CheckpointReporter(checkpoint_file, checkpoint_interval))

    # Add a progress bar for the long production run
    simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval * 10, step=True,
                                                   progress=True, remainingTime=True, speed=True,
                                                   totalSteps=production_steps + simulation.currentStep, separator='\t'))

    # Run production!
    simulation.step(production_steps)

    print(f"\n--- Simulation for {output_prefix} finished! ---")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 run_openmm_simulation.py <input_pdb_file> <output_prefix>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    prefix = sys.argv[2]

    run_simulation(pdb_file, prefix)
