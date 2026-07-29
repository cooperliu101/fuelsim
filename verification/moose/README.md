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
mesh counts, heat-source continuation, and normal penalty coefficient. The
fuel is 10.000 mm high and the cladding is 10.020 mm high, leaving 20 um of
axial projection margin in the reference geometry.

The contact normalization is essential:

```text
formulation = penalty
penalty = 1e14
normalize_penalty = true
```

MOOSE's node-face penalty is a nodal spring. Multiplication by its RZ nodal
area makes the nodal force `K*A*g`. Fuelsim now uses the same NTS
discretization: a unique fuel-node projection, current fuel half-edge
tributary area, and equal-and-opposite force distributed by the cladding line
shape functions. A MOOSE mortar-penalty input with the same numeric `penalty`
does not have this meaning and is not an equivalent oracle.

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
  Outputs/file_base=/tmp/fuelsim_m1_nts/m1 \
  Outputs/exodus=false \
  Outputs/console=false
```

The checked output was byte-for-byte identical to
`m1_fuel_cladding_gap_rz_out.csv`. The final fuel-surface vector postprocessor
output is preserved as `m1_fuel_surface_final.csv`. The scalar final row
supplies three temperatures and three displacements, while the surface file
supplies all 11 node coordinates, displacements, contact pressures, nodal
areas, and penetrations used by `tests/solver_tests.cpp`.

The final fuelsim-to-MOOSE differences are:

```text
contact pressure relative L2 error:       0.2207%
maximum nodal pressure relative error:    0.3328%
total contact force relative error:       0.0443%
projected / active fuel surface nodes:    11 / 11
```

All six temperature/displacement differences are also below 1%. Contact
pressure is therefore an acceptance metric rather than a diagnostic-only
quantity.
