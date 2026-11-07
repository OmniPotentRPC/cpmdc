#!/usr/bin/env python3
"""Independent structure-factor sum for the scaffold oracle.

The C agreement test freezes the printed energy and forces.
"""

import math


def main() -> int:
    npw = 8
    cell = 10.0
    alpha = 0.5
    sigma = 1.0
    pos = (0.2, -0.3, 0.4)
    zval = 3
    n = 1
    while n * n * n < npw:
        n += 2
    half = n // 2
    energy = 0.0
    forces = [0.0, 0.0, 0.0]
    g_visited = 0
    atom_visits = 0
    for k in range(npw):
        x = k % n
        y = (k // n) % n
        zidx = k // (n * n)
        ix = x - half
        iy = y - half
        iz = zidx - half
        gx = 2.0 * math.pi / cell * ix
        gy = 2.0 * math.pi / cell * iy
        gz = 2.0 * math.pi / cell * iz
        g2 = (gx * gx + gy * gy) + gz * gz
        g_visited += 1
        if g2 == 0.0:
            continue
        weight = (4.0 * math.pi / g2) * math.exp(-g2 / (4.0 * alpha)) * math.exp(
            -g2 / (2.0 * sigma * sigma)
        )
        phase = gx * pos[0] + gy * pos[1] + gz * pos[2]
        wz = weight * zval
        energy += wz * math.cos(phase)
        s = math.sin(phase)
        forces[0] += wz * s * gx
        forces[1] += wz * s * gy
        forces[2] += wz * s * gz
        atom_visits += 1
    print(f"{energy:.16e}")
    print(f"{forces[0]:.16e} {forces[1]:.16e} {forces[2]:.16e}")
    print(g_visited, atom_visits)
    energy_ref = 1.7274791563162685e02
    force_ref = (8.1165490211079074e00, -7.4606561327671184e00, 2.7471849298200780e01)
    if abs(energy - energy_ref) >= 1e-8:
        raise SystemExit(f"oracle energy {energy:.16e} ref {energy_ref:.16e}")
    for got, ref in zip(forces, force_ref):
        if abs(got - ref) >= 1e-8:
            raise SystemExit(f"oracle force {got:.16e} ref {ref:.16e}")
    if g_visited != 8 or atom_visits != 8:
        raise SystemExit(f"oracle visits {g_visited} {atom_visits}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
