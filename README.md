# Conditional Modal Decomposition of Aerodynamic Flow Fields

This repository contains a Python implementation of a conditional modal decomposition framework for unsteady aerodynamic flow analysis. The methodology extracts coherent flow structures associated with lift and drag dynamics using conditional averaging and modal projection techniques applied to pressure and velocity fields.The workflow is designed for CFD post-processing of large datasets generated from unsteady simulations. For more details please visit:
https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/structure-and-dynamics-of-the-laminar-separation-bubble/8DA7E999B511882B64AC059181D738FB

---

# Methodology

The instantaneous fluctuating flow field is defined as:

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

contains the instantaneous pressure and velocity components, and

$$
\Phi_{mean}
$$

represents the time-averaged flow field.

---

# Conditional Averaging and Modal Construction

The framework constructs conditional aerodynamic modes associated with the dominant lift and drag dynamics of the flow.

The lift-related conditional mode:

$$
\phi_{cl}
$$

is computed using snapshots satisfying:

$$
C_L > \overline{C_L}
$$

where

$$
\overline{C_L}
$$

is the mean lift coefficient.

Similarly, the drag-related conditional mode:

$$
\phi_{cd}
$$

is constructed from residual fluctuating fields satisfying:

$$
C_D > \overline{C_D}
$$

where

$$
\overline{C_D}
$$

is the mean drag coefficient.

This conditional averaging approach isolates coherent flow structures associated with high-lift and high-drag events during the oscillation cycle.

---

# Modal Projection

The framework computes the temporal modal coefficients:

$$
a_{cl}, \quad a_{cd}
$$

through inner-product projections of the fluctuating flow fields onto the conditional modes.

The lift modal coefficient is obtained from:

$$
a_{cl} = \langle \phi_{cl}, \Phi' \rangle
$$

After removing the contribution of the lift mode, the residual fluctuating field is defined as:

$$
\Phi_{res} =

\Phi'-

a_{cl}\phi_{cl}
$$

The drag modal coefficient is then computed from:

$$
a_{cd} = \langle \phi_{cd}, \Phi_{res} \rangle
$$

This decomposition separates the dominant lift-related dynamics from the remaining residual flow structures associated with drag fluctuations.

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
└── calculate_acd.py
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



# Input Files and Grid Configuration

The repository expects the following CFD input files:

```text
mesh_file.x
mean_flow.mat
cl.out
cd.out
snapshot-*
```

These files contain the computational mesh, mean flow solution, aerodynamic coefficients, and instantaneous flow snapshots required for the conditional modal decomposition workflow.

---

# Structured Grid Definition

The structured Pointwise mesh dimensions are defined in `main.py` as:

```python
ni = 601
nj = 301
```

where:

* `ni` is the maximum number of indices in the streamwise direction
* `nj` is the maximum number of indices in the wall-normal direction

The original CFD mesh dimensions are defined as:

```python
ni_ansys = ni - 1
nj_ansys = nj - 1
```

The CFD solution is interpolated from the original solver mesh onto the structured Pointwise mesh using:

$$
q_{ordered} = W q_{unordered}
$$

where:

* $q_{unordered}$ is the original CFD field
* $q_{ordered}$ is the structured Pointwise field
* $W$ is the sparse interpolation matrix

---

# Snapshot Synchronization and Sampling

The aerodynamic coefficients and flow snapshots are synchronized using:

```python
start_index = 24004
```

This index corresponds to the beginning of the analyzed oscillation cycle used for the modal decomposition.

The flow snapshots are sampled using:

```python
snapshot_step = 5
```

This means the flow variables were stored every 5 timesteps relative to the aerodynamic coefficient history files (`CL` and `CD`) to reduce memory and storage requirements.

The aerodynamic coefficients are sampled as:

```python
cl = cl1[start_index::snapshot_step, 1]
cd = cd1[start_index::snapshot_step, 1]
```

to ensure synchronization with the stored CFD snapshots.

---

# `mesh_file.x`

Structured Pointwise mesh file containing the computational grid coordinates.

The file stores:

* Streamwise coordinates (x)
* Wall-normal coordinates (y)

which are reshaped into the structured mesh:

$$
(x_{grid}, y_{grid})
$$

with dimensions:

$$
(ni, nj)
$$

This mesh defines the target interpolation domain used throughout the modal decomposition process.

---

# `mean_flow.mat`

MATLAB binary file containing the mean flow solution exported from the CFD solver.

The file contains the mean flow state:

$$
\Phi_{mean} =
\begin{bmatrix}
p_{mean} \
u_{mean} \
v_{mean}
\end{bmatrix}
$$

stored on the original CFD mesh.

The variables are interpolated onto the structured Pointwise mesh before computing fluctuating quantities.

---

# `cl.out`

Time history file of the lift coefficient:

$$
C_L(t)
$$

This file is used for:

* Computing the mean lift coefficient

$$
\overline{C_L}
$$

* Selecting conditional snapshots satisfying:

$$
C_L > \overline{C_L}
$$

for construction of the lift-related conditional mode:

$$
\phi_{cl}
$$

---

# `cd.out`

Time history file of the drag coefficient:

$$
C_D(t)
$$

This file is used for:

* Computing the mean drag coefficient

$$
\overline{C_D}
$$

* Selecting conditional residual snapshots satisfying:

$$
C_D > \overline{C_D}
$$

for construction of the drag-related conditional mode:

$$
\phi_{cd}
$$

---

# `snapshot-*`

Instantaneous CFD flow snapshots exported from the CFD solver.

Each snapshot contains:

$$
\Phi =
\begin{bmatrix}
p \
u \
v
\end{bmatrix}
$$

defined on the original CFD mesh.

The snapshots are automatically loaded and sorted using Python `glob`.

Each instantaneous field undergoes the following operations:

1. Interpolation onto the structured Pointwise mesh
2. Mean subtraction
3. Construction of fluctuating fields

$$
\Phi' = \Phi - \Phi_{mean}
$$

4. Projection onto the conditional modes:

$$
a_{cl} = \langle \phi_{cl}, \Phi' \rangle
$$

$$
a_{cd} = \langle \phi_{cd}, \Phi_{res} \rangle
$$

where:

$$
\Phi_{res}
$$

is the residual fluctuating field after removing the lift-related contribution.
---

# Dependencies

The code requires:

* NumPy
* SciPy
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
