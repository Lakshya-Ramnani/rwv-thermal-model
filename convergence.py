"""
Mesh convergence study for the RWV insert thermal model.

Refines the (Nr, Nz) grid at the baseline configuration and compares
transition time and uniformity gap against a fine reference solution.

The explicit timestep is recomputed for every grid as 40% of the stability
limit, so this varies dr, dz and dt together rather than dr and dz alone.

Author: Lakshya Ramnani
"""

import numpy as np
import matplotlib.pyplot as plt

from rwv_model import BASELINE, solve

BLUE, CORAL = '#185FA5', '#993C1D'

GRIDS = [(6, 9), (8, 13), (10, 15), (12, 21), (16, 29), (20, 37), (28, 53), (40, 75)]
REFERENCE = (56, 105)


def study(grids=GRIDS, reference=REFERENCE, t_max=1500.0):
    ref = solve(BASELINE, t_max=t_max, Nr=reference[0], Nz=reference[1])
    print(f"Reference grid {reference[0]}x{reference[1]}: "
          f"transition {ref['transition'] / 60:.5f} min, gap {ref['gap']:.2f} s\n")

    print(f"{'grid':>10} {'transition (min)':>18} {'gap (s)':>10} {'rel. error':>12}")
    rows = []
    for Nr, Nz in grids:
        r = solve(BASELINE, t_max=t_max, Nr=Nr, Nz=Nz)
        err = 100 * (r['transition'] - ref['transition']) / ref['transition']
        rows.append((Nr * Nz, r['transition'] / 60, r['gap'], abs(err)))
        print(f"{Nr:>4}x{Nz:<5} {r['transition'] / 60:>18.5f} "
              f"{r['gap']:>10.2f} {err:>11.3f}%")
    return np.array(rows), ref


if __name__ == '__main__':
    rows, ref = study()
    nodes, times, gaps, errs = rows.T

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    ax1.plot(nodes, times, 'o-', color=BLUE, lw=2, ms=6)
    ax1.axhline(ref['transition'] / 60, color=CORAL, ls='--', lw=1.2,
                label=f"Reference {REFERENCE[0]}x{REFERENCE[1]}")
    ax1.set_xscale('log')
    ax1.set_xlabel('Total grid nodes (Nr x Nz)', fontsize=10)
    ax1.set_ylabel('Transition time (min)', fontsize=10)
    ax1.set_title('Grid convergence: transition time', fontsize=11)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    ax2.loglog(nodes, errs, 's-', color=CORAL, lw=2, ms=6)
    ax2.set_xlabel('Total grid nodes (Nr x Nz)', fontsize=10)
    ax2.set_ylabel('Absolute relative error (%)', fontsize=10)
    ax2.set_title('Convergence of transition time vs reference', fontsize=11)
    ax2.grid(True, alpha=0.3, which='both')

    fig.suptitle('Mesh Convergence, Perforated-Channel RWV Insert', fontsize=13)
    fig.tight_layout()
    fig.savefig('Figure6_convergence.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print("\nSaved Figure6_convergence.png")
