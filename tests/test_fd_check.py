#!/usr/bin/env python3
"""Run cpmdc-fd-check on the reference evaluator.

The reference forces are the analytic gradient of its energy, so central
differences must agree. Its periodic stress is -dE/d(eps) / volume, the
OpenCPMD sign, so --stress agrees as well. The numerical virials vyy and
vzz equal -dE/d(eps) of
E = sum_i 0.5 k Z_i |r_i|^2 (Bohr) + 1e-10 sum |h|^2 (Angstrom), with
k = 1e-3 * 40 / 70 and a 10 Angstrom cube. An isolated cell has no stress.
"""

import subprocess
import sys

BOHR = 0.529177210903
WATER = [(8, (0.0, 0.0, 0.1173)), (1, (0.0, 0.7572, -0.4692)), (1, (0.0, -0.7572, -0.4692))]


def main() -> int:
    tool, isolated, periodic, step = sys.argv[1:5]
    tight = ["--max-relative-error", "1e-6"]
    cases = [
        ("forces", 0, "PASS", [*tight, isolated, step]),
        ("forces", 0, "PASS", [*tight, periodic, step]),
        ("virial", 0, "PASS", ["--stress", periodic, step]),
        ("stress not computed", 2, "", ["--stress", isolated, step]),
    ]
    failures = []
    for prefix, want_rc, verdict, args in cases:
        proc = subprocess.run([tool, *args], capture_output=True, text=True, check=False)
        out = proc.stdout + proc.stderr
        last = [ln for ln in out.splitlines() if ln.startswith(prefix)][-1:]
        if proc.returncode != want_rc or not last or not last[0].endswith(verdict):
            failures.append(f"{args}: exit {proc.returncode}, want {want_rc}\n{out}")
        elif " checked 0 of " in last[0]:
            failures.append(f"{args}: no component was checked\n{out}")
        print(last[0] if last else out)
        if prefix == "virial":
            for axis, label in ((1, "vyy "), (2, "vzz ")):
                d = sum(1e-3 * 40 / 70 * z * (p[axis] / BOHR) ** 2 for z, p in WATER) + 2e-8
                got = [float(ln.split()[1]) for ln in out.splitlines() if ln.startswith(label)]
                if len(got) != 1 or abs(got[0] + d) > 1e-9 * d:
                    failures.append(f"{label}numerical {got}, want {-d:.12e}")
    print("\n".join(failures), file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
