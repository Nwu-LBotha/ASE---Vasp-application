from ase.build import bulk
from ase.calculators.vasp import Vasp
import os

# This part searches if matplot lib is installed. 
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Note: matplotlib not found. Plots will be skipped.")

###############################
# SETUP
###############################

def get_silicon():
    """Create silicon structure."""
    return bulk('Si', 'diamond', a=5.43)


def run_calculation(label, atoms, encut, kpts):
    """Run a single VASP calculation and return energy per atom."""

    atoms = atoms.copy()

    calc = Vasp(
        directory=label,

        xc='PBE',
        encut=encut,
        kpts=kpts,

        ismear=0,      ## 0 is more suited for semiconductors
        sigma=0.05,

        lwave=False,
        lcharg=False,
    )

    atoms.calc = calc
    energy = atoms.get_potential_energy()

    return energy / len(atoms)


# Ensure output directory exists
os.makedirs("results/convergence", exist_ok=True)


###############################
# ENCUT CONVERGENCE
###############################

print("=" * 60)
print("ENCUT Convergence Test")
print("=" * 60)
print()

encut_values = [200, 250, 300, 350, 400, 450, 500]
encut_energies = []

for encut in encut_values:
    atoms = get_silicon()
    label = f"results/convergence/encut_{encut}"

    energy = run_calculation(label, atoms, encut, (4, 4, 4))
    encut_energies.append(energy)

    print(f"ENCUT = {encut:4d} eV  ->  E = {energy:.6f} eV/atom")

# Reference = highest cutoff
encut_ref = encut_energies[-1]
encut_diff = [(e - encut_ref) * 1000 for e in encut_energies]

print()
print("Energy differences (meV/atom):")
for encut, diff in zip(encut_values, encut_diff):
    print(f"ENCUT = {encut:4d} eV  ->  ΔE = {diff:+8.2f} meV/atom")

for i, diff in enumerate(encut_diff):
    if abs(diff) < 1.0:
        print(f"\nRecommended ENCUT: {encut_values[i]} eV")
        break


###############################
# K-POINT CONVERGENCE
###############################

print()
print("=" * 60)
print("K-point Convergence Test")
print("=" * 60)
print()

kpt_values = [2, 3, 4, 5, 6, 8]
kpt_energies = []

for k in kpt_values:
    atoms = get_silicon()
    label = f"results/convergence/kpts_{k}x{k}x{k}"

    energy = run_calculation(label, atoms, 400, (k, k, k))
    kpt_energies.append(energy)

    print(f"{k}x{k}x{k}  ->  E = {energy:.6f} eV/atom")

# Better reference (min energy)
kpt_ref = min(kpt_energies)
kpt_diff = [(e - kpt_ref) * 1000 for e in kpt_energies]

print()
print("Energy differences (meV/atom):")
for k, diff in zip(kpt_values, kpt_diff):
    print(f"{k}x{k}x{k}  ->  ΔE = {diff:+8.2f} meV/atom")

for i, diff in enumerate(kpt_diff):
    if abs(diff) < 1.0:
        print(f"\nRecommended k-points: {kpt_values[i]}x{kpt_values[i]}x{kpt_values[i]}")
        break


###############################
# PLOTTING
###############################

if HAS_MATPLOTLIB:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # ENCUT plot
    ax1.plot(encut_values, encut_diff, 'o-', linewidth=2)
    ax1.axhline(y=1, color='g', linestyle='--', label='+1 meV')
    ax1.axhline(y=-1, color='g', linestyle='--', label='-1 meV')
    ax1.axhline(y=0, color='k', linewidth=0.5)
    ax1.set_xlabel('ENCUT (eV)')
    ax1.set_ylabel('ΔE (meV/atom)')
    ax1.set_title('ENCUT Convergence')
    ax1.grid(True)
    ax1.legend()

    # k-point plot
    ax2.plot(kpt_values, kpt_diff, 's-', linewidth=2)
    ax2.axhline(y=1, color='g', linestyle='--', label='+1 meV')
    ax2.axhline(y=-1, color='g', linestyle='--', label='-1 meV')
    ax2.axhline(y=0, color='k', linewidth=0.5)
    ax2.set_xlabel('k-point grid (NxNxN)')
    ax2.set_ylabel('ΔE (meV/atom)')
    ax2.set_title('K-point Convergence')
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("convergence_plots.png", dpi=150)

    print("\nSaved: convergence_plots.png")


###############################
# END
###############################

print()
print("=" * 60)
print("Convergence testing complete!")
print("=" * 60)
print()
