# Conditional Modal Decomposition of Aerodynamic Flow Fields

This repository contains a Python implementation of a conditional modal decomposition framework for unsteady aerodynamic flow analysis.

The methodology extracts coherent flow structures associated with lift and drag dynamics using conditional averaging and modal projection techniques applied to pressure and velocity fields.

The workflow is designed for CFD post-processing of large datasets generated from unsteady simulations.

---

# Methodology

The instantaneous fluctuating flow field is defined as

$$
\Phi' = \Phi - \Phi_{mean}
$$

where

$$
\Phi =
\begin{bmatrix}
p \
u \
v
\end{bmatrix}
$$

contains pressure and velocity components.

The framework computes:

Lift-related conditional mode:

$$
\phi_{cl}
$$

Drag-related conditional mode:

$$
\phi_{cd}
$$

together with their temporal modal coefficients:

$$
a_{cl}, \quad a_{cd}
$$

The modal coefficients are obtained through projection:

$$
a_{cl} = \langle \phi_{cl}, \Phi' \rangle
$$

$$
a_{cd} = \langle \phi_{cd}, \Phi_{res} \rangle
$$

where $\Phi_{res}$ is the residual fluctuating field after removing the contribution of the lift mode.

---

# Main Features

* Conditional averaging based on aerodynamic coefficients
* Sparse interpolation operators for high performance
* Fast structured/unstructured mesh mapping
* Reusable interpolation matrices
* Nearest-neighbor fallback interpolation
* CFD-oriented data processing pipeline
* Memory-efficient implementation

---

# Repository Structure

```text
.
├── main.py
├── load_mean_flow.py
├── normalize_mode.py
├── interpolation.py
├── calculate_phi_cl.py
├── calculate_acl_and_phicd.py
├── calculate_phi_cd.py
├── calculate_acd.py
├── aoa_18.x
├── mean_puv_aoa18.mat
├── cl18.out
├── cd18.out
└── p18-*
```

---

# Workflow

The main routine performs the following operations:

1. Load structured mesh and mean flow
2. Build sparse interpolation operators
3. Interpolate ANSYS flow fields onto Pointwise mesh
4. Compute fluctuating fields
5. Construct lift conditional mode `phi_cl`
6. Compute lift modal coefficients `acl`
7. Construct residual flow fields
8. Build drag conditional mode `phi_cd`
9. Compute drag modal coefficients `acd`

---

# Interpolation Strategy

The interpolation is written as

[
q_{ordered} = W q_{unordered}
]

where:

* ( q_{unordered} ) represents the original ANSYS field
* ( q_{ordered} ) represents the structured Pointwise field
* ( W ) is the sparse interpolation matrix

The interpolation operator is constructed only once and reused throughout the simulation, significantly reducing computational cost.

---

# Conditional Averaging

The lift mode is computed using snapshots satisfying:

[
CL > \overline{CL}
]

Similarly, the drag mode is constructed from residual fields satisfying:

[
CD > \overline{CD}
]

This approach isolates coherent structures associated with high-lift and high-drag events.

---

# Grid Definitions

The structured mesh dimensions are defined as:

```python
ni = 601
nj = 301
```

where:

* `ni` is the maximum number of indices in the streamwise direction
* `nj` is the maximum number of indices in the wall-normal direction

The original ANSYS mesh dimensions are:

```python
ni_ansys = ni - 1
nj_ansys = nj - 1
```

---

# Snapshot Sampling

The flow-field snapshots are sampled using:

```python
snapshot_step = 5
```

This means flow variables were stored every 5 timesteps compared to the aerodynamic coefficient files (`CL` and `CD`) to reduce memory requirements.

The synchronization starting location is defined using:

```python
start_index = 24004
```

which corresponds to the beginning of the analyzed oscillation cycle.

---

# Dependencies

The code requires:

* NumPy
* SciPy

Install using:

```bash
pip install numpy scipy
```

---

# Running the Code

Execute the main script:

```bash
python main.py
```

---

# Output

The code computes:

* `phi_cl`
* `phi_cd`
* `acl`
* `acd`

together with diagnostic information such as:

* minimum/maximum modal coefficients
* NaN detection
* number of conditional samples

Optional saving of modes and coefficients can be added using:

```python
np.save(...)
```

---

# Applications

This framework can be used for:

* Low-Reynolds-number aerodynamics
* Laminar separation bubble analysis
* Dynamic stall studies
* Modal decomposition
* Coherent structure identification
* Reduced-order modeling
* CFD post-processing

---

# Notes

* Physics and interpolation behavior are preserved relative to the original implementation.
* Sparse interpolation matrices are reused to improve performance.
* Nearest-neighbor interpolation is only applied outside the Delaunay triangulation region.

---
