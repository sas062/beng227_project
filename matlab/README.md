# Ren et al. 2022 α-β phase-locking model (MATLAB)

This folder contains a MATLAB implementation of a coupled phase-oscillator
model for pancreatic α and β cells, following the mathematical structure used
in:

Ren, H., Li, Y., Han, C. *et al.* Pancreatic α and β cells are globally
phase-locked. **Nat Commun** 13, 3721 (2022).

## Files

- `ren2022_default_params.m` – baseline parameters and simulation setup.
- `ren2022_rhs.m` – ODE right-hand side implementing coupled phase equations.
- `run_ren2022_demo.m` – runs the simulation, plots trajectories, and reports
  an estimated winding-number ratio.

## Usage

From MATLAB:

```matlab
cd matlab
out = run_ren2022_demo();
```

Override parameters (example):

```matlab
p = struct('K_alpha_to_beta', 0.35, 'K_beta_to_alpha', -0.55, ...
           'tspan', [0, 5400]);
out = run_ren2022_demo(p);
```

## Notes

- The code is designed so coupling signs match the interpretation from the
  paper: α→β is stimulatory and β→α is inhibitory.
- If you have experimentally fitted parameters or exact supplementary values,
  replace fields in `ren2022_default_params.m` or pass overrides to
  `run_ren2022_demo`.
