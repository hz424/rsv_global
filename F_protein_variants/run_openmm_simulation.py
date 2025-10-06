#!/usr/bin/env python3
"""
Run an MD simulation of a protein complex with OpenMM.
Fixed for antibody–antigen work: consistent FF/water, 0.15 M salt, CUDA props, resumable.
"""

import sys, os
from pathlib import Path
from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer

def _ensure_parent_dirs(*paths):
    """Create parent directories for all given paths if they don't exist."""
    for p in paths:
        d = Path(p).parent
        if str(d) not in ("", "."):
            d.mkdir(parents=True, exist_ok=True)

def run_simulation(pdb_filename, output_prefix):
    print("="*80)
    print(f"  OpenMM MD Simulation")
    print("="*80)
    print(f"Input PDB:       {pdb_filename}")
    print(f"Output prefix:   {output_prefix}")
    print(f"Process ID:      {os.getpid()}")
    print(f"SLURM_PROCID:    {os.environ.get('SLURM_PROCID', 'N/A')}")
    print(f"SLURM_LOCALID:   {os.environ.get('SLURM_LOCALID', 'N/A')}")
    print("="*80)

    # --- Simulation Parameters ---
    temperature = 300 * kelvin
    time_step   = 0.004 * picoseconds   # 4 fs with HMR
    equilibration_steps = 125_000       # 0.5 ns
    production_steps    = 12_500_000    # 50 ns (extended from 20 ns)
    report_every        = 2_500
    checkpoint_every    = 12_500
    salt_molar          = 0.15 * molar

    eq_ns  = (equilibration_steps * time_step).value_in_unit(nanoseconds)
    prod_ns= (production_steps    * time_step).value_in_unit(nanoseconds)
    print(f"\nSimulation Parameters:")
    print(f"  Temperature:     {temperature}")
    print(f"  Time step:       {time_step.value_in_unit(femtoseconds):.0f} fs (with HMR)")
    print(f"  Equilibration:   {eq_ns:.2f} ns ({equilibration_steps:,} steps)")
    print(f"  Production:      {prod_ns:.0f} ns ({production_steps:,} steps)")
    print(f"  Salt:            {salt_molar}")
    print(f"  Report interval: every {report_every} steps")

    # --- Output paths ---
    traj_path   = f"{output_prefix}_trajectory.dcd"
    cpt_path    = f"{output_prefix}_checkpoint.chk"
    min_pdb     = f"{output_prefix}_minimized.pdb"
    init_solv_pdb = f"{output_prefix}_initial_solvated.pdb"
    system_xml  = f"{output_prefix}_system.xml"
    integrator_xml = f"{output_prefix}_integrator.xml"
    log_path    = f"{output_prefix}_log.txt"

    # Ensure all parent directories exist
    _ensure_parent_dirs(traj_path, cpt_path, min_pdb, init_solv_pdb,
                        system_xml, integrator_xml, log_path)

    print(f"\nOutput files:")
    print(f"  Trajectory:      {traj_path}")
    print(f"  Checkpoint:      {cpt_path}")
    print(f"  Minimized PDB:   {min_pdb}")
    print(f"  Log:             {log_path}")

    # If a trajectory already exists, we'll append to it on resume.

    # --- 1) Build / Load System ---
    print("\n1) Building system...")
    # Force field + water model (standard TIP3P for compatibility)
    forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')

    # Check if we can resume from saved system XML (avoids water bond issues)
    if os.path.exists(system_xml):
        print(f"   > Resuming from saved system: {system_xml}")
        with open(system_xml, 'r') as f:
            system = XmlSerializer.deserialize(f.read())
        # Load topology from initial solvated PDB
        pdb = PDBFile(init_solv_pdb)
        topology, positions = pdb.topology, pdb.positions
    elif os.path.exists(cpt_path) and os.path.exists(init_solv_pdb):
        # Checkpoint exists AND init_solv_pdb exists - recreate system from it
        # (init_solv_pdb may have water bond issues, but we'll save system XML to avoid this in future)
        print(f"   > Checkpoint exists, recreating system from initial solvated PDB")
        pdb = PDBFile(init_solv_pdb)
        topology, positions = pdb.topology, pdb.positions
        print("   > Creating System (PME, HMR, HBond constraints, barostat)...")
        system = forcefield.createSystem(topology,
                                         nonbondedMethod=PME,
                                         nonbondedCutoff=1.0*nanometers,
                                         constraints=HBonds,
                                         hydrogenMass=4*amu)    # HMR
        # Save system XML for future restarts
        with open(system_xml, 'w') as f:
            f.write(XmlSerializer.serialize(system))
        print(f"   > Saved system to: {system_xml}")
    else:
        # Fresh start - prepare system from scratch
        print("   > No checkpoint. Preparing from input PDB...")
        fixer = PDBFixer(filename=pdb_filename)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        # IMPORTANT: do NOT strip glycans by default; uncomment next line only if you intend to remove ALL HETATMs.
        # fixer.removeHeterogens(True)   # keepWater=True, but removes ligands/glycans too
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)

        modeller = Modeller(fixer.topology, fixer.positions)
        print("   > Adding solvent and 0.15 M salt...")
        modeller.addSolvent(forcefield, padding=1.2*nanometers, model='tip3p',
                            ionicStrength=salt_molar)
        topology  = modeller.topology
        positions = modeller.positions
        with open(init_solv_pdb, 'w') as fh:
            PDBFile.writeFile(topology, positions, fh)

        print("   > Creating System (PME, HMR, HBond constraints, barostat)...")
        system = forcefield.createSystem(topology,
                                         nonbondedMethod=PME,
                                         nonbondedCutoff=1.0*nanometers,
                                         constraints=HBonds,
                                         hydrogenMass=4*amu)    # HMR
        # Save system XML for future restarts (avoids water bond issues)
        with open(system_xml, 'w') as f:
            f.write(XmlSerializer.serialize(system))
        print(f"   > Saved system to: {system_xml}")

    system.addForce(MonteCarloBarostat(1*bar, temperature, 25))

    # --- 2) Platform / Integrator ---
    print("\n2) Selecting platform...")
    integrator = LangevinMiddleIntegrator(temperature, 1/picosecond, time_step)

    # Stable, reproducible RNG seed
    seed_env = (os.environ.get('SLURM_ARRAY_TASK_ID') or
                os.environ.get('SLURM_PROCID') or "0")
    seed = int(seed_env) + 12345
    try:
        integrator.setRandomNumberSeed(seed)
    except Exception:
        pass

    properties = {'CudaPrecision': 'mixed'}
    gpu_env = os.environ.get('CUDA_VISIBLE_DEVICES', None)
    if gpu_env:
        # When CUDA_VISIBLE_DEVICES is set, always use device index 0
        # (the first device in the visible list)
        properties['CudaDeviceIndex'] = '0'
        print(f"   > GPU Assignment: Using first visible GPU (CUDA_VISIBLE_DEVICES = {gpu_env})")
        print(f"   > OpenMM will use device index 0")
    else:
        print("   > WARNING: CUDA_VISIBLE_DEVICES not set; using OpenMM default GPU selection")

    try:
        platform   = Platform.getPlatformByName('CUDA')
        simulation = Simulation(topology, system, integrator, platform, properties)
        print(f"   > ✓ Successfully initialized CUDA platform")
        print(f"   > Platform: {platform.getName()}")
        print(f"   > Precision: {properties.get('CudaPrecision', 'default')}")
    except Exception as e:
        print(f"   > ✗ CUDA unavailable ({e})")
        print(f"   > Falling back to CPU platform (will be SLOW!)")
        platform   = Platform.getPlatformByName('CPU')
        simulation = Simulation(topology, system, integrator, platform)

    simulation.context.setPositions(positions)

    # --- 3) Minimize / 4) Equilibrate (if fresh) ---
    if not os.path.exists(cpt_path):
        print("\n3) Minimizing...")
        simulation.minimizeEnergy()
        min_positions = simulation.context.getState(getPositions=True).getPositions()
        with open(min_pdb, 'w') as fh:
            PDBFile.writeFile(simulation.topology, min_positions, fh)
        print(f"   > Saved minimized coordinates: {min_pdb}")

        print(f"\n4) Equilibrating for {eq_ns:.2f} ns...")
        simulation.reporters.append(StateDataReporter(sys.stdout, report_every,
                                step=True, potentialEnergy=True, temperature=True,
                                progress=True, remainingTime=True, speed=True,
                                totalSteps=equilibration_steps, separator='\t'))
        simulation.step(equilibration_steps)
        simulation.saveCheckpoint(cpt_path)
        print(f"   > Equilibration done. Checkpoint: {cpt_path}")

    # --- 5) Production ---
    print(f"\n5) Production: {prod_ns:.0f} ns...")
    if os.path.exists(cpt_path):
        simulation.loadCheckpoint(cpt_path)
        # Derive progress from saved time (currentStep isn't restored by checkpoints)
        curr_ps = simulation.context.getState().getTime().value_in_unit(picoseconds)
        dt_ps = time_step.value_in_unit(picoseconds)
        current_step = int(round(curr_ps / dt_ps))
        print(f"   > Loaded checkpoint at t = {curr_ps:.3f} ps (step {current_step})")
    else:
        current_step = 0

    # Calculate remaining steps
    steps_completed = current_step - equilibration_steps if current_step > equilibration_steps else 0
    steps_remaining = production_steps - steps_completed

    if steps_remaining <= 0:
        print(f"   > Production already complete! ({steps_completed} steps done)")
        print("\n" + "="*80)
        print(f"  SIMULATION ALREADY COMPLETED")
        print("="*80)
        return

    print(f"   > Steps completed: {steps_completed}/{production_steps}")
    print(f"   > Steps remaining: {steps_remaining}")

    # reporters
    # Determine if we're appending (resuming) or starting fresh
    is_resuming = (current_step > 0)

    # Safer: only append to files that actually exist
    dcd_append = is_resuming and os.path.exists(traj_path)
    log_append = is_resuming and os.path.exists(log_path)

    simulation.reporters.clear()
    simulation.reporters.append(DCDReporter(traj_path, report_every, append=dcd_append))
    simulation.reporters.append(StateDataReporter(log_path, report_every,
                                step=True, time=True, potentialEnergy=True,
                                temperature=True, volume=True, speed=True, append=log_append))
    simulation.reporters.append(CheckpointReporter(cpt_path, checkpoint_every))
    simulation.reporters.append(StateDataReporter(sys.stdout, report_every*10,
                                step=True, progress=True, remainingTime=True,
                                speed=True, totalSteps=steps_remaining,
                                separator='\t'))

    simulation.step(steps_remaining)

    print("\n" + "="*80)
    print(f"  SIMULATION COMPLETED SUCCESSFULLY")
    print("="*80)
    print(f"Output prefix:   {output_prefix}")
    print(f"Trajectory:      {traj_path}")
    print(f"Log file:        {log_path}")
    print(f"Checkpoint:      {cpt_path}")
    print(f"Total time:      {eq_ns:.2f} ns (equil) + {prod_ns:.0f} ns (prod) = {eq_ns + prod_ns:.2f} ns")
    print("="*80)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_openmm_simulation.py <input_pdb_file> <output_prefix>")
        sys.exit(1)
    run_simulation(sys.argv[1], sys.argv[2])
