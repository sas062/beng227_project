# Hybrid Alpha-Beta Islet Model

This project implements a hybrid ODE-PDE model of pancreatic islet dynamics in MATLAB. Alpha and beta cells are represented as phase oscillators, while glucagon is modeled as a spatially distributed field on a 2D grid.

The model extends a reduced alpha-beta phase framework by adding:
- Spatial organization of alpha and beta cells.
- Heterogeneous intrinsic beta-cell frequencies.
- A glucagon diffusion-reaction PDE.
- Spatially weighted beta-to-alpha coupling.
- Local glucagon sensing by beta cells through bilinear interpolation.

## Model overview

The system contains three coupled components:

1. Alpha-cell phase dynamics.
2. Beta-cell phase dynamics.
3. A 2D glucagon field \(G(x,y,t)\).

### Alpha-cell dynamics

Each alpha cell evolves according to an intrinsic frequency plus beta-derived input:

\[
\frac{d\theta_\alpha}{dt}
=
\omega_\alpha
+
K_{\beta\alpha}\, I_\beta \, f_{r\alpha}(\theta_\alpha)
\]

where the effective beta-cell signal received by each alpha cell is a weighted spatial average of beta-cell secretion:

\[
I_{\beta,n}
=
\sum_{m=1}^{N_\beta}
W_{nm}\, f_s(\theta_m^\beta)
\]

with distance-based weights

\[
\widetilde{W}_{nm}
=
\exp\!\left(
-\frac{\|x_n^\alpha - x_m^\beta\|^2}{2\ell_I^2}
\right),
\qquad
W_{nm}
=
\frac{\widetilde{W}_{nm}}{\sum_{k=1}^{N_\beta}\widetilde{W}_{nk}}
\]

### Beta-cell dynamics

Each beta cell evolves according to its intrinsic frequency plus glucagon-dependent coupling:

\[
\frac{d\theta_\beta}{dt}
=
\omega_\beta
+
K_{\alpha\beta}\, G(x_m^\beta,t)\, f_{r\beta}(\theta_\beta)
\]

The glucagon value at each beta-cell location is obtained from the grid by bilinear interpolation.

### Glucagon PDE

The glucagon field evolves as:

\[
\frac{dG}{dt}
=
D_G \nabla^2 G
+
S(x,y,t)
-
\tau G
\]

where the secretion source is

\[
S(x,y,t)
=
\sum_{n=1}^{N_\alpha}
s_G\, f_s(\theta_n^\alpha(t))\, \phi(x-x_n^\alpha, y-y_n^\alpha)
\]

and \(\phi\) is a Gaussian kernel centered at each alpha-cell location.

## Spatial organization

The computational domain is a 2D square grid. Cells are placed in a circular islet:
- Beta cells are biased```markdown
# Hybrid Alpha-Beta Islet Model

This project implements a hybrid ODE-PDE model of pancreatic islet dynamics in MATLAB. Alpha and beta cells are represented as phase oscillators, while glucagon is modeled as a spatially distributed field on a 2D grid.

The model extends a reduced alpha-beta phase framework by adding:
- Spatial organization of alpha and beta cells.
- Heterogeneous intrinsic beta-cell frequencies.
- A glucagon diffusion-reaction PDE.
- Spatially weighted beta-to-alpha coupling.
- Local glucagon sensing by beta cells through bilinear interpolation.

## Model Overview

The system contains three coupled components:

1. Alpha-cell phase dynamics.
2. Beta-cell phase dynamics.
3. A 2D glucagon field \(G(x,y,t)\).

### Alpha-cell dynamics

Each alpha cell evolves according to an intrinsic frequency plus beta-derived input:

\[
\frac{d\theta_\alpha}{dt}
=
\omega_\alpha
+
K_{\beta\alpha}\, I_\beta \, f_{r\alpha}(\theta_\alpha)
\]

where the effective beta-cell signal received by each alpha cell is a weighted spatial average of beta-cell secretion:

\[
I_{\beta,n}
=
\sum_{m=1}^{N_\beta}
W_{nm}\, f_s(\theta_m^\beta)
\]

with distance-based weights

\[
\widetilde{W}_{nm}
=
\exp\!\left(
-\frac{\|x_n^\alpha - x_m^\beta\|^2}{2\ell_I^2}
\right),
\qquad
W_{nm}
=
\frac{\widetilde{W}_{nm}}{\sum_{k=1}^{N_\beta}\widetilde{W}_{nk}}
\]

### Beta-cell dynamics

Each beta cell evolves according to its intrinsic frequency plus glucagon-dependent coupling:

\[
\frac{d\theta_\beta}{dt}
=
\omega_\beta
+
K_{\alpha\beta}\, G(x_m^\beta,t)\, f_{r\beta}(\theta_\beta)
\]

The glucagon value at each beta-cell location is obtained from the grid by bilinear interpolation.

### Glucagon PDE

The glucagon field evolves as:

\[
\frac{dG}{dt}
=
D_G \nabla^2 G
+
S(x,y,t)
-
\tau G
\]

where the secretion source is

\[
S(x,y,t)
=
\sum_{n=1}^{N_\alpha}
s_G\, f_s(\theta_n^\alpha(t))\, \phi(x-x_n^\alpha, y-y_n^\alpha)
\]

and \(\phi\) is a Gaussian kernel centered at each alpha-cell location.

## Spatial Organization

The computational domain is a 2D square grid. Cells are placed in a circular islet:
- Beta cells are biased toward the interior.
- Alpha cells are biased toward the periphery.

This creates a more realistic islet architecture than uniform random placement.

## Files

### `main.m`

Sets parameters, cell counts, initial phases, grid geometry, islet layout, initial glucagon field, initial state vector, and calls the ODE solver.

### `rhs_islet.m`

Main right-hand-side function passed to the solver. It:
- Unpacks alpha phases, beta phases, and glucagon field.
- Computes beta-to-alpha signaling.
- Interpolates glucagon to beta-cell positions.
- Computes alpha and beta phase derivatives.
- Calls the glucagon PDE helper.
- Reassembles the full state derivative vector.

### `glucagon_pde_rhs.m`

Computes the PDE contribution for the glucagon field:
- Diffusion.
- Secretion source.
- Linear decay.

### `glucagon_source_term.m`

Builds the spatial glucagon secretion field from alpha-cell phases and locations using Gaussian kernels.

### `insulin_secretion.m`

Computes the effective beta-cell secretion signal received by alpha cells using spatially weighted beta-to-alpha coupling.

## Parameter Structure

The model uses a parameter struct `p` containing quantities such as:
- `p.w_a`: alpha-cell intrinsic frequency.
- `p.w_b`: heterogeneous beta-cell intrinsic frequencies.
- `p.K_ab`: glucagon-to-beta coupling coefficient.
- `p.K_ba`: beta-to-alpha coupling coefficient.
- `p.tau`: glucagon clearance rate.
- `p.DG`: glucagon diffusivity.
- `p.sG`: secretion scaling factor.
- `p.sigmaG`: Gaussian source width.
- `p.lI`: interaction length scale for beta-to-alpha weighting.
- `p.fs`, `p.fra`, `p.frb`: secretion/response functions.
- Grid and cell-position data.

## Initial Conditions

The initial state consists of:
- Random alpha phases.
- Random beta phases.
- Zero glucagon everywhere on the grid.

These are stacked into a single solver state vector:

```matlab
y0 = [theta_a0; theta_b0; G0(:)];
