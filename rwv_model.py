"""
RWV insert thermal model.

Ti-6Al-4V insert with perforated radial channels holding Pluronic F-127
thermosensitive plugs, each retained by a sintered porous titanium frit.

Physics
    Ti wall    lumped capacitance (Bi = 0.071), analytical exponential
    F-127 plug 2D axisymmetric heat equation, explicit finite differences

Plug boundary conditions
    r = R_ch   Dirichlet, T = T_Ti(t)          (k_Ti / k_f = 17.5)
    r = 0      symmetry, dT/dr = 0
    z = 0      Robin, convective loss to the cold drug chamber
    z = L_f    Robin, convective gain from warm medium through the frit

Final configuration: Ti = 5 mm, plug = 2 mm x R = 1 mm, T_transition = 31 C.
Predicted transition 4.33 min, drug release 6.3 to 7.3 min.

Author: Lakshya Ramnani
"""

import numpy as np

# Material properties: Ti-6Al-4V, Pluronic F-127
RHO_TI, CP_TI, K_TI = 4430.0, 526.0, 7.0
RHO_F, CP_F, K_F = 1050.0, 3800.0, 0.4
ALPHA_F = K_F / (RHO_F * CP_F)

BASELINE = dict(
    L_ti=5.0e-3, L_f=2.0e-3, R_ch=1.0e-3,
    T_medium=37.0, T_initial=4.0, T_drug=4.0, T_transition=31.0,
    h_medium=100.0, h_drug=10.0, h_med_face=10.0,
)


def wall(p):
    """Lumped Ti wall steady state and time constant."""
    h_tot = p['h_medium'] + p['h_drug']
    T_inf = (p['h_medium'] * p['T_medium'] + p['h_drug'] * p['T_drug']) / h_tot
    return T_inf, RHO_TI * CP_TI * p['L_ti'] / h_tot


def _laplacian(T, dr, dz, inv_r, p):
    """
    Axisymmetric Laplacian, radial part plus axial part.

    inv_r holds 1/r for the interior rows only, shaped for broadcasting.
    The last radial row is Dirichlet, so its value here is never used.
    """
    lap_r = np.zeros_like(T)
    lap_z = np.zeros_like(T)

    # Radial. At r = 0 the singular 1/r term resolves to a second d2T/dr2.
    lap_r[1:-1] = ((T[2:] - 2 * T[1:-1] + T[:-2]) / dr ** 2
                   + (T[2:] - T[:-2]) / (2 * dr) * inv_r)
    lap_r[0] = 4 * (T[1] - T[0]) / dr ** 2

    # Axial. Face terms come from eliminating the ghost node in the flux balance.
    lap_z[:, 1:-1] = (T[:, 2:] - 2 * T[:, 1:-1] + T[:, :-2]) / dz ** 2
    lap_z[:, 0] = (2 * (T[:, 1] - T[:, 0]) / dz ** 2
                   - 2 * p['h_drug'] * (T[:, 0] - p['T_drug']) / (K_F * dz))
    lap_z[:, -1] = (2 * (T[:, -2] - T[:, -1]) / dz ** 2
                    - 2 * p['h_med_face'] * (T[:, -1] - p['T_medium']) / (K_F * dz))

    return lap_r + lap_z


def solve(p=BASELINE, t_max=600.0, Nr=12, Nz=21, snap_times=(), record=False):
    """
    March the plug forward in time until every node exceeds T_transition.

    Returns first_crossing (hottest node crosses), transition (coldest node
    crosses, the true opening) and their difference as the uniformity gap.
    Set record=True to also return centreline histories and 2D snapshots.
    """
    dr, dz = p['R_ch'] / (Nr - 1), p['L_f'] / (Nz - 1)
    dt = 0.4 * 0.5 / (ALPHA_F * (1 / dr ** 2 + 1 / dz ** 2))   # 40% of the explicit limit
    T_inf, tau = wall(p)

    if T_inf < p['T_transition']:                    # wall never gets hot enough
        return dict(T_inf=T_inf, tau=tau, dr=dr, dz=dz, dt=dt, feasible=False,
                    transition=None, first_crossing=None, gap=None)

    inv_r = 1 / np.linspace(0, p['R_ch'], Nr)[1:-1, None]
    T = np.full((Nr, Nz), p['T_initial'])

    keys = ('t', 'T_ti', 'drug_face', 'centre', 'med_face')
    hist = []
    snaps = []
    pending = list(snap_times)
    record_every = max(1, int(0.5 / dt))             # steps per half second
    first = transition = None

    for step in range(int(t_max / dt) + 1):
        t = step * dt
        T[-1] = T_inf + (p['T_initial'] - T_inf) * np.exp(-t / tau)

        if first is None and T.max() >= p['T_transition']:
            first = t
        if transition is None and T.min() >= p['T_transition']:
            transition = t
            if not record:
                break

        if record:
            if step % record_every == 0:
                hist.append((t, T[-1, 0], T[0, 0], T[0, Nz // 2], T[0, -1]))
            while pending and t >= pending[0]:
                snaps.append((pending.pop(0), T.copy()))

        # The wall row is Dirichlet, so it is excluded from the update.
        T[:-1] += dt * ALPHA_F * _laplacian(T, dr, dz, inv_r, p)[:-1]

    out = dict(T_inf=T_inf, tau=tau, dr=dr, dz=dz, dt=dt,
               transition=transition, first_crossing=first,
               feasible=transition is not None,
               gap=None if transition is None else transition - first)
    if record:
        out.update(dict(zip(keys, np.array(hist).T)), snaps=snaps)
    return out


# ============================================================================
# FIGURES: final configuration
# ============================================================================
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    SNAPS = (60.0, 180.0, 260.0, 420.0)
    LABELS = tuple(f't = {s / 60:.1f} min' for s in SNAPS)
    COLORS = ('#185FA5', '#BA7517', '#3B6D11', '#A32D2D')
    Nr, Nz = 12, 21
    p = BASELINE
    Tt = p['T_transition']

    r = solve(p, t_max=600.0, Nr=Nr, Nz=Nz, snap_times=SNAPS, record=True)

    print(f"Ti steady state      {r['T_inf']:.2f} C")
    print(f"Ti time constant     {r['tau']:.1f} s ({r['tau'] / 60:.2f} min)")
    print(f"First point crosses  {r['first_crossing']:.1f} s ({r['first_crossing'] / 60:.2f} min)")
    print(f"Barrier transitions  {r['transition']:.1f} s ({r['transition'] / 60:.2f} min)")
    print(f"Uniformity gap       {r['gap']:.1f} s")
    print(f"Drug-side final T    {r['drug_face'][-1]:.2f} C")
    print(f"Grid                 dr = {r['dr'] * 1e3:.4f} mm, dz = {r['dz'] * 1e3:.4f} mm, "
          f"dt = {r['dt'] * 1e3:.3f} ms")
    print(f"Bi (Ti, outer)       {p['h_medium'] * p['L_ti'] / K_TI:.3f}  (lumped if < 0.1)")
    print(f"k_Ti / k_F-127       {K_TI / K_F:.1f}  (justifies the Dirichlet wall)")

    def finish(ax, xlabel, ylabel, title, loc='lower right'):
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=11)
        ax.legend(loc=loc, fontsize=9, framealpha=0.95)
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('white')

    geom = (f"Ti = {p['L_ti'] * 1e3:.0f} mm, plug = {p['L_f'] * 1e3:.0f} mm "
            f"x R = {p['R_ch'] * 1e3:.1f} mm, T$_{{transition}}$ = {Tt:.0f}$\\degree$C")

    # Figure 1: thermal evolution
    curves = [
        ('T_ti',      '#185FA5', 'Ti wall (lumped)'),
        ('med_face',  '#993C1D', 'F-127 medium-side face (r=0)'),
        ('centre',    '#3B6D11', 'F-127 axial center (r=0)'),
        ('drug_face', '#993556', 'F-127 drug-side face (r=0, slowest)'),
    ]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for key, color, lab in curves:
        ax.plot(r['t'] / 60, r[key], color=color, linewidth=2, label=lab)
    ax.axvline(5.0, color='#444441', ls='--', lw=1.2, alpha=0.8, label='Spin-up complete (5 min)')
    ax.axhline(Tt, color='#A32D2D', ls='--', lw=1.2, alpha=0.7,
               label=f'F-127 transition ({Tt:.0f}$\\degree$C)')
    ax.axvline(r['transition'] / 60, color='#3C3489', ls=':', lw=1.5,
               label=f"Barrier transitions ({r['transition'] / 60:.2f} min)")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 38)
    finish(ax, 'Time (minutes)', 'Temperature ($\\degree$C)',
           f'Thermal evolution: final configuration\n{geom}')
    fig.savefig('figure_thermal_evolution.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    # Figure 2: 2D temperature field
    fig, axes = plt.subplots(1, 4, figsize=(15, 4), constrained_layout=True)
    for ax, (_, Ts), lab in zip(axes, r['snaps'], LABELS):
        im = ax.imshow(Ts.T, origin='lower', aspect='auto', cmap='RdBu_r', vmin=4, vmax=37,
                       extent=[0, p['R_ch'] * 1e3, 0, p['L_f'] * 1e3])
        ax.contour(np.linspace(0, p['R_ch'] * 1e3, Nr), np.linspace(0, p['L_f'] * 1e3, Nz),
                   Ts.T, levels=[Tt], colors='black', linewidths=1.5)
        ax.set_title(lab, fontsize=11)
        ax.set_xlabel('Radial position r (mm)', fontsize=10)
    axes[0].set_ylabel('Axial position z (mm)\n$\\leftarrow$ drug face       medium face $\\rightarrow$',
                       fontsize=10)
    fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02, label='Temperature ($\\degree$C)')
    fig.suptitle(f'2D temperature field in F-127 plug. '
                 f'Black contour: T = {Tt:.0f}$\\degree$C transition isotherm', fontsize=11, y=1.02)
    fig.savefig('figure_2d_field.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    # Figure 3: axial centreline profile
    fig, ax = plt.subplots(figsize=(10, 5.5))
    z_mm = np.linspace(0, p['L_f'] * 1e3, Nz)
    for (t_s, Ts), color in zip(r['snaps'], COLORS):
        ax.plot(z_mm, Ts[0], color=color, linewidth=2, marker='o', markersize=4,
                label=f't = {t_s / 60:.1f} min')
    ax.axhline(Tt, color='#A32D2D', ls='--', lw=1.2, alpha=0.7,
               label=f'F-127 transition ({Tt:.0f}$\\degree$C)')
    ax.set_xlim(0, p['L_f'] * 1e3)
    ax.set_ylim(0, 38)
    finish(ax, 'Axial position in F-127 plug (mm)\n'
               '$\\leftarrow$ drug-chamber face          medium-side face $\\rightarrow$',
           'Temperature at centerline r = 0 ($\\degree$C)',
           f'Axial temperature profile: F-127 plug centerline\n{geom}')
    fig.savefig('figure_axial_profile.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    print("\nSaved figure_thermal_evolution.png, figure_2d_field.png, figure_axial_profile.png")
