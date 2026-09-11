"""
Parameter sensitivity for the RWV insert thermal model.

Design knobs, set by hardware and recipe
    T_transition  F-127 concentration
    L_ti          Ti wall thickness
    L_f           F-127 plug axial length
    R_ch          channel radius

Uncertainty parameters, set by the RWV and the operating protocol
    h_medium      outer surface in gently rotating medium
    h_med_face    tissue-side plug face, behind the frit
    h_drug        drug-side plug face and inner Ti surface

Sweeps run at the final 31 C baseline so the figures match the reported
configuration. Coarse grid (Nr=10, Nz=15) is within 0.04% of a 
56x105 reference and is sufficient for trends.

Author: Lakshya Ramnani
"""

import numpy as np
import matplotlib.pyplot as plt

from rwv_model import BASELINE, solve

BLUE, CORAL = '#185FA5', '#993C1D'
COARSE = dict(Nr=10, Nz=15)      # sweep grid, within 0.04% of a 56x105 reference
FINE = dict(Nr=12, Nz=21)        # reporting grid, matches rwv_model.py


def sweep(key, values, label, t_max=1500.0):
    """Vary one parameter, return transition time (min) and uniformity gap (s)."""
    print(f"\n{label}")
    times, gaps = [], []
    for v in values:
        r = solve({**BASELINE, key: v}, t_max=t_max, **COARSE)
        times.append(r['transition'] / 60 if r['feasible'] else np.nan)
        gaps.append(r['gap'] if r['feasible'] else np.nan)
        print(f"  {v:>9.4g}  transition {times[-1]:6.2f} min   gap {gaps[-1]:5.1f} s")
    return np.array(times), np.array(gaps)


def floor_span(ax, values, min_span):
    """Stop a near-flat result from being magnified into an apparent trend."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size and np.ptp(finite) < min_span:
        mid = 0.5 * (finite.min() + finite.max())
        ax.set_ylim(mid - min_span / 2, mid + min_span / 2)


def panel(ax, x, times, gaps, xlabel, title, baseline):
    ax2 = ax.twinx()
    ax.plot(x, times, 'o-', color=BLUE, lw=2, ms=6)
    ax2.plot(x, gaps, 's--', color=CORAL, lw=1.5, ms=5)
    ax.axvline(baseline, color='gray', ls=':', alpha=0.6)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel('Transition time (min)', color=BLUE, fontsize=10)
    ax2.set_ylabel('Uniformity gap (s)', color=CORAL, fontsize=10)
    ax.tick_params(axis='y', labelcolor=BLUE)
    ax2.tick_params(axis='y', labelcolor=CORAL)
    ax.set_title(title, fontsize=11)
    ax.grid(True, alpha=0.3)
    floor_span(ax, times, 0.5)       # minutes
    floor_span(ax2, gaps, 5.0)       # seconds


if __name__ == '__main__':
    b = solve(BASELINE, **FINE)
    print(f"Baseline: transition {b['transition'] / 60:.2f} min, "
          f"gap {b['gap']:.1f} s, T_Ti_inf {b['T_inf']:.2f} C, tau {b['tau']:.1f} s")

    # --- Design knobs -------------------------------------------------------
    knobs = [
        ('T_transition', np.array([23, 25, 27, 29, 31, 33]),
         'Transition temperature ($\\degree$C)',
         'Knob 1: F-127 concentration (transition temperature)', 31, 1.0),
        ('L_ti', np.array([2, 3, 4, 5, 6, 7, 8, 10]) * 1e-3,
         'Ti wall thickness (mm)', 'Knob 2: Ti wall thickness', 5.0, 1e3),
        ('L_f', np.array([1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]) * 1e-3,
         'F-127 plug axial length (mm)', 'Knob 3: F-127 plug length', 2.0, 1e3),
        ('R_ch', np.array([0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5]) * 1e-3,
         'Channel radius (mm)', 'Knob 4: Channel radius', 1.0, 1e3),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for ax, (key, vals, xlabel, title, base, scale) in zip(axes.flat, knobs):
        times, gaps = sweep(key, vals, title)
        panel(ax, vals * scale, times, gaps, xlabel, title, base)
    fig.suptitle('Design Knob Sensitivity, Perforated-Channel RWV Insert', fontsize=13, y=1.00)
    fig.tight_layout()
    fig.savefig('figureA_design_knobs.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    # --- Uncertainty parameters --------------------------------------------
    uncert = [
        ('h_medium', np.array([50, 75, 100, 150, 200, 300]),
         'h$_{medium}$ (W/m$^2\\cdot$K)', 'Outer convection (RWV medium)', 100),
        ('h_med_face', np.array([5, 10, 20, 50, 100]),
         'h$_{med,face}$ (W/m$^2\\cdot$K)', 'Tissue-side face (through frit)', 10),
        ('h_drug', np.array([2, 5, 10, 14, 16, 18]),
         'h$_{drug}$ (W/m$^2\\cdot$K)', 'Drug-side (stagnant chamber)', 10),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (key, vals, xlabel, title, base) in zip(axes, uncert):
        times, gaps = sweep(key, vals, title, t_max=3000.0)
        panel(ax, vals, times, gaps, xlabel, title, base)
    fig.suptitle('Uncertainty Parameter Sensitivity, Convective Coefficients', fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig('figureB_uncertainty.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    # --- Robustness and the h_drug failure threshold ------------------------
    print("\nRobustness, h_medium +/- 25%")
    for f in (0.75, 1.00, 1.25):
        h = BASELINE['h_medium'] * f
        r = solve({**BASELINE, 'h_medium': h}, t_max=1500.0)
        print(f"  h_medium = {h:5.1f}  transition {r['transition'] / 60:.2f} min")

    # Wall steady state equals T_transition at the critical h_drug.
    # Solving (h_m T_m + h_d T_d) / (h_m + h_d) = T_t for h_d:
    p = BASELINE
    h_crit = (p['h_medium'] * (p['T_medium'] - p['T_transition'])
              / (p['T_transition'] - p['T_drug']))
    print(f"\nCritical h_drug (device never gels above this): {h_crit:.1f} W/m^2 K")

    print("\nSaved figureA_design_knobs.png, figureB_uncertainty.png")
