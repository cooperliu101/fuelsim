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

## M2.1 transient heat capacity

`m21_transient_heat_rz.i` isolates the consistent transient heat-capacity
residual on a 4x2 Quad4 RZ mesh. The cylinder is insulated, has uniform
properties, and is heated by a constant volumetric source:

```text
density = 10000 kg/m^3
specific_heat = 300 J/(kg K)
heat_source = 3e6 W/m^3
initial_temperature = 600 K
```

The exact spatially uniform solution is `T(t) = 600 + t` K. The MOOSE objects
used by this reference are:

```text
ADHeatConduction
ADHeatConductionTimeDerivative
ADBodyForce
ADGenericConstantMaterial
```

Both thermal kernels explicitly use the reference mesh. The time integrator is
implicit Euler with `dt=1 s` and `end_time=10 s`. Reproduce the checked run
from this directory with:

```bash
source /home/cooper/miniforge/etc/profile.d/conda.sh
conda activate moose
/home/cooper/projects/july/july-opt \
  -i m21_transient_heat_rz.i \
  Outputs/file_base=/tmp/fuelsim_m21_transient_heat_rz \
  Outputs/console=false
```

`m21_transient_heat_rz_out.csv` preserves the complete 11-row temperature
history. At the final time:

```text
MOOSE average temperature:        610 K
analytic average temperature:     610 K
relative error:                   0%
MOOSE nodal L2 error:             0
acceptance threshold:             < 0.1%
```

This case verifies the heat-capacity term and time integration. A nonuniform
manufactured solution is still required before claiming spatial transient
conduction verification.

## M2.2 Norton creep

`m22_norton_creep_rz.i` is a one-Quad4 homogeneous RZ material-point proxy.
It applies a constant 100 MPa axial tensile traction to an isotropic cylinder
with `E=200 GPa` and `nu=0.3`. The generic secondary Norton law uses:

```text
coefficient = 1e-30 Pa^-3 s^-1
n_exponent = 3
m_exponent = 0
activation_energy = 0
```

The MOOSE constitutive chain is:

```text
ADPowerLawCreepStressUpdate
  -> ADComputeMultipleInelasticStress
  -> AD small incremental strain mechanics
```

`ADMaterialRealAux` exposes the stateful `effective_creep_strain` property.
Using the non-AD `MaterialRealAux` object is not equivalent and produces an
AD/non-AD material-property error. Reproduce the checked run with:

```bash
source /home/cooper/miniforge/etc/profile.d/conda.sh
conda activate moose
/home/cooper/projects/july/july-opt \
  -i m22_norton_creep_rz.i \
  Outputs/file_base=/tmp/fuelsim_m22_norton_creep_rz \
  Outputs/console=false
```

For the nominal constant stress, the analytic values at `t=100 s` are
`effective_creep_strain=1e-4` and
`axial_displacement=0.001*(100e6/200e9 + 1e-4)=6e-7 m`.
The final row of `m22_norton_creep_rz_out.csv` gives:

```text
axial stress:                     99.998007620195 MPa
stress error from 100 MPa:        0.0019924%
effective creep strain:           9.9991036503199e-5
effective creep strain error:     0.0089635%
axial displacement:               5.9997808603447e-7 m
axial displacement error:         0.0036523%
acceptance threshold for each:    < 0.1%
```

The axial, radial, and hoop creep strains have the expected J2 ratio
`1 : -0.5 : -0.5`, so the creep strain is deviatoric to output precision.

## M2.2 J2 plasticity

`m22_j2_plastic_rz.i` uses the same one-Quad4 RZ material-point geometry and
applies a monotonic axial strain of 0.002. The material constants are:

```text
E = 200 GPa
nu = 0.3
yield_stress = 200 MPa
hardening_constant = 2 GPa
```

The MOOSE constitutive chain is:

```text
ADIsotropicPlasticityStressUpdate
  -> ADComputeMultipleInelasticStress
  -> AD small incremental strain mechanics
```

Reproduce the checked run with:

```bash
source /home/cooper/miniforge/etc/profile.d/conda.sh
conda activate moose
/home/cooper/projects/july/july-opt \
  -i m22_j2_plastic_rz.i \
  Outputs/file_base=/tmp/fuelsim_m22_j2_plastic_rz \
  Outputs/console=false
```

For linear isotropic hardening, the one-dimensional analytic final state is:

```text
effective plastic strain:
  (E*strain - yield_stress)/(E + hardening_constant)
  = 9.9009900990099e-4

axial stress:
  E*(strain - effective_plastic_strain)
  = 201.9801980198 MPa
```

The final row of `m22_j2_plastic_rz_out.csv` matches both values to the printed
precision:

```text
MOOSE effective plastic strain:   9.9009900990099e-4
MOOSE axial stress:               201.9801980198 MPa
relative errors:                  < 1e-12%
acceptance threshold for each:    < 0.1%
```

The axial, radial, and hoop plastic strains also have the J2 ratio
`1 : -0.5 : -0.5`. This monotonic case does not replace a future unload/reload
history test.

## M2 reference environment and conventions

All three M2 reference inputs were syntax-checked and solved with one MPI rank
and one thread using:

```text
MOOSE commit:          93b11698be3fcd33049ae73e32f411fb2985261d
July commit:           a96d73792bee7c5f54eb65e33b04487b24276a27
Executable SHA256:     1cb3a0fbf5650addd087ceeb8521f82d2ddb11650c0932b7ec274226430e0ee4
Result:                every time step reported Solve Converged!
```

The July worktree was dirty, so the executable hash remains part of the
provenance. For these RZ cases MOOSE component names map as follows:

```text
xx = radial
yy = axial
zz = hoop
xy = rz shear
```

Both inelastic models require incremental strain and committed quadrature-point
history. Trial history must be recomputed from the same old state during every
global Newton or line-search evaluation and committed only after a converged
time step. If creep and plasticity are later combined in one MOOSE reference,
the required model order is `inelastic_models = 'creep plasticity'`.
