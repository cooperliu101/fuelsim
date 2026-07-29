# MOOSE references

## M0 single fuel cylinder

`m0_simple_fuel_rz.i` is the independent reference for the M0 steady
thermoelastic fuel cylinder. Both heat conduction and mechanics use the
reference mesh so that the weak form matches fuelsim exactly.

The checked run used:

```text
July/MOOSE executable: /home/cooper/projects/july/july-opt
MOOSE commit:          93b11698be
PETSc:                 3.25.2
MPI ranks / threads:   1 / 1
Mesh:                  40 radial x 10 axial Quad4
DOFs:                  1353
Result:                Solve Converged!
```

Reproduce from this directory:

```bash
source /home/cooper/miniforge/etc/profile.d/conda.sh
conda activate moose
/home/cooper/projects/july/july-opt -i m0_simple_fuel_rz.i
```

The final row of `m0_simple_fuel_rz_out.csv` is the acceptance snapshot used by
`tests/solver_tests.cpp`.

## M1 fuel, cladding, gap heat, and contact

`m1_fuel_cladding_gap_rz.i` matches the default fuelsim M1 geometry, materials,
mesh counts, heat-source continuation, and normal penalty coefficient.

The contact normalization is essential:

```text
formulation = penalty
penalty = 1e14
normalize_penalty = true
```

MOOSE's node-face penalty is a nodal spring. Multiplication by its RZ nodal
area makes the nodal force `K*A*g`, which is the lumped counterpart of
fuelsim's interface traction `p=K*max(-g,0)`. A MOOSE mortar-penalty input with
the same numeric `penalty` does not have this meaning and is not an equivalent
oracle.

The checked run used:

```text
July/MOOSE executable: /home/cooper/projects/july/july-opt
MOOSE commit:          93b11698be3fcd33049ae73e32f411fb2985261d
July commit:           a96d73792bee7c5f54eb65e33b04487b24276a27
Executable SHA256:     1cb3a0fbf5650addd087ceeb8521f82d2ddb11650c0932b7ec274226430e0ee4
PETSc / SLEPc:         3.25.2 / 3.25.0
MPI ranks / threads:   1 / 1
Mesh:                  fuel 40x10, cladding 6x10 Quad4
Nodes / elements:      528 / 460
Nonlinear DOFs:        1584
Loading:               20 steps, dt=0.05, heat source proportional to time
Result:                all 20 steps reported Solve Converged!
```

The July worktree was dirty, so the executable hash is recorded in addition to
the commit. The exact checked command was:

```bash
source /home/cooper/miniforge/etc/profile.d/conda.sh
conda activate moose
/home/cooper/projects/july/july-opt \
  -i m1_fuel_cladding_gap_rz.i \
  Outputs/file_base=/tmp/fuelsim_m1_final_check \
  Outputs/exodus=false
```

The checked output was byte-for-byte identical to
`m1_fuel_cladding_gap_rz_out.csv`. Its final row supplies six M1 acceptance
values in `tests/solver_tests.cpp`: three temperatures and three
displacements. All six fuelsim differences are below 1%.

Peak node pressure is diagnostic rather than an acceptance scalar. On this
mesh MOOSE reports 5.450439 MPa while fuelsim reports about 6.63 MPa; the
different endpoint/contact interpolation makes the peak mesh-sensitive.
Future pressure validation should use a refined-mesh total load or
area-weighted norm.
