# Ren et al. 2022 α-β phase-locking model (MATLAB)

This repository contains a MATLAB implementation of a coupled phase-oscillator
model for pancreatic α and β cells, following the mathematical structure used
in:

Ren, H., Li, Y., Han, C. *et al.* Pancreatic α and β cells are globally
phase-locked. **Nat Commun** 13, 3721 (2022).

## Files (repo root)

- `ren2022_default_params.m` – baseline parameters and simulation setup.
- `ren2022_rhs.m` – ODE right-hand side implementing coupled phase equations.
- `run_ren2022_demo.m` – runs the simulation, plots trajectories, and reports
  an estimated winding-number ratio.
- `run_ren2022_paperstyle_figures.m` – renders paper-style validation panels
  (activity traces, relative phase, locking-angle distribution, event raster,
  and a locking-mode map).

## Usage

From MATLAB (opened at repository root):

```matlab
out = run_ren2022_demo();
```

Override parameters (example):

```matlab
p = struct('K_alpha_to_beta', 0.35, 'K_beta_to_alpha', -0.55, ...
           'tspan', [0, 5400]);
out = run_ren2022_demo(p);
```

## Why these plots?

The demo figures are chosen as **diagnostic plots for phase locking**, not as a
pixel-match recreation of a specific figure panel from the paper:

1. **Wrapped phases** (`φα`, `φβ` vs time): checks both oscillators are active
   and shows timing relationship.
2. **Relative phase** (`φα - φβ` vs time): most direct readout of phase locking;
   bounded/near-constant behavior indicates locking.
3. **Relative phase histogram**: reveals the preferred phase offset (locking
   angle) over the full simulation.
4. **Cycle counts + winding ratio**: verifies global locking ratio via the
   empirical `n_alpha / n_beta` winding number estimate.

## Paper-style visual validation

To generate figure panels designed for quick visual comparison against the
paper's phase-locking style of analysis, run:

```matlab
out = run_ren2022_paperstyle_figures();
```

This creates:
- time traces (activity proxies for α/β),
- relative phase trajectory and polar histogram,
- cycle-event raster from phase wraps, and
- a period-sweep locking-mode map with current parameters marked.

## Notes

- The code is designed so coupling signs match the interpretation from the
  paper: α→β is stimulatory and β→α is inhibitory.
- If you have experimentally fitted parameters or exact supplementary values,
  replace fields in `ren2022_default_params.m` or pass overrides to
  `run_ren2022_demo`.
