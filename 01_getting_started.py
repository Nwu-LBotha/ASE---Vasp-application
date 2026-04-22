#Please not that this python code should be in the running folder while the job is submitted in the job script.
# The jobscript should the one used to run jobs on either your computer or on a super computer.  
# You will need a license for VASP

from ase.build import bulk
from ase.calculators.vasp import Vasp

# Build structure
atoms = bulk('Si', 'diamond', a=5.43)

print("Structure created:")
print(f"  Formula: {atoms.get_chemical_formula()}")
print(f"  Number of atoms: {len(atoms)}")
print(f"  Lattice constant: {atoms.cell[0, 1]:.3f} Å")
print()

# Set up calculator
calc = Vasp(
    directory='01_si_energy_results',

    xc='PBE',
    encut=300,
    kpts=(4, 4, 4),

    ismear=1,
    sigma=0.1,

    lwave=False,
    lcharg=False,
)

# Attach calculator
atoms.calc = calc

print("Running VASP calculation...")
print("=" * 50)

# Run calculation
energy = atoms.get_potential_energy()

print("=" * 50)
print()

print("Results:")
print(f"  Total energy: {energy:.6f} eV")
print(f"  Energy per atom: {energy / len(atoms):.6f} eV/atom")

# Fermi level (safe way)
try:
    fermi = calc.get_fermi_level()
    print(f"  Fermi level: {fermi:.4f} eV")
except:
    print("  Fermi level: Not available")

print()
print(f"Output files written to: {calc.directory}/")

# Access additional results
fermi = calc.results.get('fermi_level')
if fermi:
    print(f"  Fermi level: {fermi:.4f} eV")

# Check convergence
if calc.results.get('converged', True):
    print("  Status: Converged successfully")
else:
    print("  Status: WARNING - Not converged!")

print()
print(f"Output files written to: {calc.directory}/")
print()
print("Congratulations! Your first VASP calculation is complete.")
print("Next: Try 02_convergence/ to learn about parameter convergence.")
