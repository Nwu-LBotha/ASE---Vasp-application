import numpy as np
from ase.build import bulk
from ase.calculators.vasp import Vasp

print("=" * 60)
print("Structure Relaxation Example")
print("=" * 60)
print()

##############################################################################
# Create distorted silicon structure
##############################################################################

# Start with experimental lattice constant (we'll let VASP find the PBE value)
atoms = bulk('Si', 'diamond', a=5.43)

# Slightly distort the structure to give the optimizer something to do
# Move one atom off its ideal position
atoms.positions[1] += [0.1, 0.05, -0.08]
atoms.pbc = True  # This is required for VASP

print("Initial structure:")
print(f"  Lattice constant: {atoms.cell[0, 1]:.4f} Å")
print(f"  Atom 1 position: {atoms.positions[0]}")
print(f"  Atom 2 position: {atoms.positions[1]}")
print()

##############################################################################
# Part 1: Ionic Relaxation (ISIF = 2)
##############################################################################

print("=" * 60)
print("Part 1: Ionic Relaxation (ISIF=2)")
print("=" * 60)
print()

# Copy structure (IMPORTANT: separate run)
atoms_ionic = atoms.copy()

calc_ionic = Vasp(
    directory='relax_ionic',
    xc='PBE',
    encut=400,
    kpts=(4, 4, 4),
    ismear=1,
    sigma=0.1,

    ibrion=2,
    isif=2,
    nsw=50,
    ediffg=-0.02,

    lwave=False,
    lcharg=False,
)

atoms_ionic.calc = calc_ionic

energy_ionic = atoms_ionic.get_potential_energy()
forces_ionic = atoms_ionic.get_forces()
max_force_ionic = np.max(np.abs(forces_ionic))

print("Results after ionic relaxation:")
print(f"  Energy: {energy_ionic:.6f} eV")
print(f"  Maximum force: {max_force_ionic:.4f} eV/Å")
print("  Relaxed positions:")
print(f"    Atom 1: {atoms_ionic.positions[0]}")
print(f"    Atom 2: {atoms_ionic.positions[1]}")
print()

##############################################################################
# Part 2: Full Cell Relaxation (ISIF = 3)
##############################################################################

print("=" * 60)
print("Part 2: Full Cell Relaxation (ISIF=3)")
print("=" * 60)
print()

# Start from slightly wrong lattice constant
atoms_full = bulk('Si', 'diamond', a=5.50)
atoms_full.pbc = True  # REQUIRED

calc_full = Vasp(
    directory='relax_full',
    xc='PBE',
    encut=400,
    kpts=(4, 4, 4),
    ismear=1,
    sigma=0.1,

    ibrion=2,
    isif=3,
    nsw=50,
    ediffg=-0.02,

    lwave=False,
    lcharg=False,
)

atoms_full.calc = calc_full

energy_full = atoms_full.get_potential_energy()
forces_full = atoms_full.get_forces()
max_force_full = np.max(np.abs(forces_full))

cell = atoms_full.get_cell()

print("Results after full relaxation:")
print(f"  Energy: {energy_full:.6f} eV")
print(f"  Maximum force: {max_force_full:.4f} eV/Å")
print(f"  Relaxed cell:\n{cell}")
print()


##############################################################################
# Summary
##############################################################################

print("=" * 60)
print("Summary")
print("=" * 60)
print()
print("Relaxation type comparison:")
print(f"  Ionic only (ISIF=2): E = {energy_ionic:.4f} eV")
print(f"  Full cell (ISIF=3):  E = {energy_full:.4f} eV")
print()
print("Key points:")
print("  - Use ISIF=2 for surfaces (where the cell size is fixed)")
print("  - Use ISIF=3 for bulk structures")
print("  - PBE typically overestimates lattice constants by ~1%")
print("  - Check OUTCAR to verify convergence")
print()
print("Next: Try 04_equation_of_state/ to calculate bulk modulus.")
