# Hybrid Alpha-Beta Islet Model (MATLAB)

This repository contains MATLAB simulations for pancreatic islet dynamics, combining:
- phase-oscillator models for alpha and beta cells, and
- a spatial glucagon field governed by a diffusion-reaction PDE.

The codebase includes both:
1. a **Recreation of the Initial Model by Ren et al.** (`ren_model.m`) used to reproduce the classic alpha-beta phase behavior, and
2. a **Hybrid Spatial ODE-PDE Model** (`main.m` + helper functions) with multiple cells on a 2D grid.

## Repository Contents

- `main.m`  
  Entry point for the spatial hybrid model. Defines parameters, initializes cell positions and phases, builds the glucagon grid state, and runs `ode15s`.

- `rhs_islet.m`  
  Right-hand-side function for the spatial model ODE system. Unpacks the state vector into alpha phases, beta phases, and glucagon field; computes coupled dynamics; repacks derivatives.

- `insulin_secretion.m`  
  Computes beta-to-alpha signaling using distance-weighted coupling between beta and alpha cell coordinates.

- `glucagon_secretion.m`  
  Constructs the spatial glucagon source term as a sum of Gaussian kernels centered at alpha-cell locations.

- `glucagon_pde_rhs.m`  
  Computes glucagon PDE terms (diffusion + secretion - linear decay).

- `ren_model.m`  
  Reduced non-spatial reference model (single alpha, single beta, scalar glucagon) that generates coupling-regime phase plots.

- `get_NN`
  Creates connectivity matrix for beta-beta interactions.

- `assign_cell_type`
  Assigns cell type to the islet structure.

- `Figures/`  
  Saved figures (PNG/FIG) for representative coupling regimes.

## Mathematical Structure

### 1) Alpha-cell phase dynamics

For each alpha cell:

$\frac{d\theta_\alpha}{dt}=\omega_\alpha + K_{\beta\alpha}\,I_\beta\,f_{r\alpha}(\theta_\alpha)$

where \(I_\beta\) is a weighted spatial average of beta secretion.

### 2) Beta-cell phase dynamics

For each beta cell:

$\frac{d\theta_\beta}{dt}=\omega_\beta + K_{\alpha\beta}\,G(x_\beta,y_\beta,t)\,f_{r\beta}(\theta_\beta)$

with glucagon sampled from the grid via bilinear interpolation.

### 3) Glucagon PDE

$\frac{\partial G}{\partial t}=D_G\nabla^2G + S(x,y,t) - \tau G$

where \(S\) is the alpha-cell secretion field built from normalized Gaussian kernels.

## Default model setup in `main.m`

- **Simulation window:** `tspan = [0 360]`
- **Cell counts:** `Na = 15`, `Nb = 85`
- **Domain:** `Lx = Ly = 400` µm
- **Grid:** `Nx = Ny = 81`
- **Architecture:** beta cells biased toward islet interior, alpha cells toward periphery
- **Solver:** `ode15s` with `RelTol=1e-6`, `AbsTol=1e-8`

## How to run

### Spatial hybrid model

From MATLAB, run:

```matlab
main
```
```matlab
plot_results_V2
```

This computes the coupled alpha/beta/glucagon dynamics and produces relevant plots to visualize dynamics.

### Reduced reference model

Run:

```matlab
ren_model
```

This generates a 4-panel figure comparing coupling regimes (`Slow`, `Fast`, `Mixed`, `Weakly Coupled`).

## Notes

- Additional development files are included in this repository. Only the files listed above in 'Repository Contents' are necessary to run the model.