# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib>=3.8,<4"]
# ///
"""Draw energy.png and scf-steps.png from scan.csv, and print the numbers
that docs/orgmode/tutorials/bond-scan.org quotes."""

import csv
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HARTREE_EV = 27.211386245988
here = pathlib.Path(__file__).resolve().parent
rows = list(csv.DictReader((here / "scan.csv").open()))
col = {k: [float(row[k]) for row in rows] for k in rows[0]}
r, energy, force = col["r_oh_ang"], col["energy_h"], col["f_bond_h_bohr"]
warm, cold = col["scf_steps"], col["scf_steps_cold"]
low = energy.index(min(energy))

for name, ylabel, series in [
    ("energy.png", "Energy above the lowest point (eV)",
     [("o-", None, [(e - energy[low]) * HARTREE_EV for e in energy])]),
    ("scf-steps.png", "SCF iterations",
     [("s-", "deck's guess (CPMDC_WARM_RMSD=0)", cold),
      ("o-", "previous call's orbitals (default)", warm)]),
]:
    fig, ax = plt.subplots(figsize=(5.0, 3.4), dpi=150)
    for style, label, y in series:
        ax.plot(r, y, style, label=label)
    ax.set_xlabel("O-H bond length (Angstrom)")
    ax.set_ylabel(ylabel)
    if name == "scf-steps.png":
        ax.set_ylim(0, max(cold) + 2)
        ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(here / name, metadata={"Software": None})

print(f"lowest r_oh={r[low]:.3f} energy_h={energy[low]:.8f}")
print(f"span_ev={(max(energy) - energy[low]) * HARTREE_EV:.3f}")
for a, b, fa, fb in zip(r, r[1:], force, force[1:]):
    if fa > 0 >= fb:
        print(f"force sign change between {a:.3f} and {b:.3f}")
print(f"scf_steps warm total={sum(warm):.0f} after_first={sorted(set(warm[1:]))}")
print(f"scf_steps cold total={sum(cold):.0f} range={min(cold):.0f}..{max(cold):.0f}")
diff = max(abs(a - b) for a, b in zip(energy, col["energy_h_cold"]))
print(f"max |energy warm - cold| = {diff:.1e} Hartree")
