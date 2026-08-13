# Solver benchmarks

`fuelsim_m1_single_core_benchmark` is the manual M1 linear-solver benchmark.
The default invocation retains the original medium case and direct solver. The
accepted arguments are:

```text
fuelsim_m1_single_core_benchmark \
  [medium|large] \
  [direct|block_jacobi|field_split|hypre] \
  [load_steps] \
  [unscaled|scaled]
```

The two meshes are:

```text
medium fuel/cladding: 100/16 radial x 64 axial
medium nodes/elements/DOFs: 7,670 / 7,424 / 23,010
large fuel/cladding:  200/32 radial x 64 axial
large nodes/elements/DOFs: 15,210 / 14,848 / 45,630
default load steps: 20
```

It is intentionally not a CTest because one run takes tens of seconds. Build
Release and run it on one pinned core:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/tmp/adlite-fuelsim-install
cmake --build build --parallel

env \
  OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  taskset -c 0 \
  ./build/fuelsim_m1_single_core_benchmark medium direct 20 unscaled
```

The matching MOOSE command is:

```bash
env \
  OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  taskset -c 0 \
  /home/cooper/projects/july/july-opt \
    -i verification/moose/m1_fuel_cladding_gap_rz.i \
    Mesh/fuel_mesh/nx=100 \
    Mesh/fuel_mesh/ny=64 \
    Mesh/clad_mesh/nx=16 \
    Mesh/clad_mesh/ny=64 \
    Outputs/csv=false \
    Outputs/exodus=false \
    Outputs/console=false
```

## 2026-07-31 measurements

Both executables were Release builds pinned to CPU 0. MOOSE used one MPI
process and one thread.

```text
default 1,584-DOF case, five-run medians:
  fuelsim optimized: 2.00 s
  MOOSE:             3.21 s
  ratio:             fuelsim 1.61x faster

default case, paired three-run medians:
  fuelsim 319e173 baseline: 2.23 s
  fuelsim optimized:        1.95 s
  object-reuse improvement: 12.6% less wall time

medium 23,010-DOF case, two-run means:
  fuelsim: 40.45 s, observed range 38.45-42.44 s
  MOOSE:   54.97 s, observed range 50.69-59.24 s
  ratio:   fuelsim 1.36x faster
```

The medium fuelsim run created one PETSc workspace for all 20 steps and
reported 62 Jacobian evaluations and 82 residual evaluations. Its warmed
internal load-path time was 37.96 s, including 14.79 s in Jacobian callbacks
and 6.12 s in residual callbacks.

These are end-to-end timings for the current reference model and direct
solver. Re-run them after changes to hardware, mesh, PETSc, linear solver,
material models, or contact algorithms.

## 2026-08-02 M3.4 measurements

All OpenMP/OpenBLAS/MKL/NumExpr thread counts were one. The 1-rank runs were
pinned to CPU 0; the 2-rank run was pinned to CPUs 0 and 1.

```text
default 1,584 DOF, 20 steps, three-run paired medians:
  697762f baseline: 0.893109 s
  M3.4:             0.872064 s
  change:           2.36% faster

medium 23,010 DOF, 20 steps:
  697762f 1 rank: 27.6641 s
  M3.4     1 rank: 27.6630 s
  M3.4     2 ranks: 10.3542 s
  observed 2-rank speedup over current 1-rank: 2.67x
```

Every medium run completed 62 nonlinear iterations with 82 residual and 62
Jacobian callbacks and created one PETSc workspace. The 2-rank result uses
distributed contribution assembly and PETSc MUMPS; it is not a claim of
general strong scaling beyond this two-rank measurement.

## 2026-08-03 M4.0 numerical-foundation check

Release builds used the same PETSc/Exodus toolchain, CPU 0, one MPI rank and
one thread for every listed numerical library. The default 1,584-DOF case was
paired against `a750379` after adding field residual diagnostics. Process-run
times were:

```text
                         first run   following two median   all-three median
  a750379 baseline:       0.874010 s       0.878382 s          0.877441 s
  M4.0 candidate:         0.882035 s       0.884105 s          0.882035 s
  all-three change:                                              0.52% slower
```

The constant-property path retains cached Lamé parameters, while nonzero
temperature coefficients use the AD-active path. The current 23,010-DOF,
20-step case completed once in `27.7143 s`, with 62 nonlinear iterations,
82 residual callbacks, 62 Jacobian callbacks and one PETSc workspace. This is
a regression check for the current benchmark only, not a broader scaling
claim.

## 2026-08-06 M5.5 shadow-state collection

The Release build used the same 23,010-DOF, 20-step case. All
OpenMP/OpenBLAS/MKL/NumExpr thread counts were one. The one-rank run was pinned
to CPU 0 and the two-rank run to CPUs 0 and 1. First-run and warmed internal
load-path times were kept separate:

```text
                         first run     warmed repeat/median
  1 rank:                27.9611 s          27.9525 s
  2 ranks:               10.5659 s          10.4864 s
```

Every run completed 62 nonlinear iterations, 82 residual callbacks, 62
Jacobian callbacks, and one PETSc workspace. These timings are current-run
observations, not a paired claim against the previous implementation and not a
general scaling result.

The old callback path stored and scattered all 23,010 state values on each
rank. With two ranks, the new contribution-derived shadow sets report:

```text
global state DOFs:                         23,010
maximum shadow DOFs on one rank:           12,290
sum of shadow DOFs across two ranks:       23,600
sum for two replicated full states:        46,020
remote shadow DOFs per callback:            8,088
remote DOFs for the old all-gather:        23,010
```

Thus the stored state-value slots summed over both ranks decrease by `48.72%`,
and the remote value payload per residual or Jacobian callback decreases from
`184,080 bytes` to `64,704 bytes`, a `64.85%` reduction. The largest rank's
persistent callback workspace contains two double buffers plus 32-bit global
indices: `245,800 bytes` instead of two replicated full double buffers totaling
`368,160 bytes`, a `33.24%` reduction. These exact counts exclude allocator,
PETSc scatter metadata, matrices, factorizations, replicated mesh geometry,
committed state, and material history; they are not process resident-set-size
claims. A single full-state gather remains after each nonlinear solve so the
existing committed-state transaction can continue on every rank, but it is no
longer performed for every callback.

## 2026-08-13 constant-time shadow-state lookup

Residual and Jacobian callbacks previously revalidated the sorted shadow
degree-of-freedom list when constructing every state view, and every contact
candidate state lookup performed a binary search. The current layout validates
the list once during solver setup and builds one contiguous 32-bit
global-to-shadow index per global degree of freedom. Callback lookups are then
direct array accesses. Missing shadow values still raise an explicit error.

The Release measurements used one process pinned to CPU 0 with OpenMP,
OpenBLAS, MKL, and NumExpr fixed to one thread. Commit `24ea1ea` is the paired
baseline. The default 1,584-degree-of-freedom case used three runs; the medium
case used one run on each revision:

```text
                                      24ea1ea       candidate      observed change
default three-run median:            1.139552 s     1.114839 s       2.17% lower
medium 23,010 DOF, 20 steps:        62.093170 s    60.722259 s       2.21% lower
```

Both revisions completed 64 iterations in the default case and 62 iterations
in the medium case, retained the same final residuals, and created one PETSc
workspace. These are paired regression measurements for the named meshes and
machine, not a general performance or scaling claim.

The lookup table costs `4 * global_state_dofs` bytes per process. The benchmark
now reports this separately as `global_to_shadow_lookup_bytes` and includes it
in `maximum_shadow_workspace_bytes`. For the 23,010-degree-of-freedom case the
additional table is 92,040 bytes. This explicit replicated index cost buys
constant-time access while the collected floating-point state values and their
per-callback communication remain limited to each process's shadow set.

## 2026-08-06 M5.6 engineering-scale iterative solvers

The Release build used two MPI processes pinned to CPUs 0 and 1. OpenMP,
OpenBLAS, MKL, and NumExpr thread counts were all one. The medium and large
cases used the same physics, 20 load steps, nonlinear tolerances, and linear
relative tolerance of `1e-8`. The direct runs used PETSc MUMPS. The iterative
runs used GMRES with the selected preconditioner. First-run and warmed internal
load-path times were recorded separately:

| DOFs | solver | completed | nonlinear / linear iterations | first | warmed |
| ---: | --- | --- | ---: | ---: | ---: |
| 23,010 | direct MUMPS | 20/20 | 62 / 62 | 10.1111 s | 10.0455 s |
| 23,010 | HYPRE | 20/20 | 84 / 7,010 | 88.0817 s | 87.6334 s |
| 45,630 | direct MUMPS | 20/20 | 62 / 62 | 21.2375 s | 21.3292 s |
| 45,630 | HYPRE | 20/20 | 77 / 2,960 | 85.9042 s | 84.4195 s |

The final global residual norms were between `2.85e-9` and `3.57e-9`, and the
existing global and per-field residual audit accepted every completed run. On
the warmed measurement, HYPRE took `8.72` times the direct time at 23,010 DOFs
and `3.96` times at 45,630 DOFs. The HYPRE iteration count decreased when the
radial resolution doubled, so these two points do not establish a monotonic
iteration or scaling trend.

The other combinations were screened on the medium mesh with two requested
load steps and the same limit of 500 linear iterations per nonlinear solve.
Neither completed the first load step:

| preconditioner | field scaling | accumulated linear iterations | final residual | internal time |
| --- | --- | ---: | ---: | ---: |
| block Jacobi | off | 1,551 | 0.7590 | 1.4133 s |
| field split | off | 2,600 | 0.6848 | 2.7615 s |
| block Jacobi | on | 1,270 | 0.9057 | 1.1073 s |
| field split | on | 1,927 | 0.9058 | 2.0678 s |

Automatic field residual scaling reduced a two-step HYPRE screen from 557 to
435 linear iterations and from 8.22 to 6.70 seconds. It did not survive the
full path: the scaled 23,010-DOF HYPRE run stopped after 13 of 20 steps with
5,185 accumulated linear iterations and a final scaled residual of `1.41`.
Consequently, the short-screen improvement is not a supported long-path
setting.

The required default 1,584-DOF direct-solver pairing used commit `9238355` as
the baseline and the same single-process CPU-0 command for the candidate. The
first internal times were `0.868126/0.878489 s`. The following two runs were
`0.869982/0.880232 s` for the baseline and `0.867705/0.881343 s` for the
candidate, giving warmed two-run medians of `0.875107/0.874524 s`. The
candidate was `0.067%` faster on that warmed statistic, much smaller than the
run-to-run spread. Both versions used 62 nonlinear and 62 linear iterations,
ended at the same `7.683024548830e-9` residual, and reported the same contact
force, so this is evidence of no material direct-path regression rather than a
speedup claim.

The matching one-process MOOSE run used the same 23,010-DOF mesh, physics,
20 steps, direct solver, CPU 0, and disabled CSV, Exodus, and console output.
Its first and warmed wall times were `47.02/46.13 s`. A warmed fuelsim run
measured with the same `/usr/bin/time` wall-clock command was `27.04 s`, which
is `41.4%` less time, or a MOOSE-to-fuelsim ratio of `1.71`. This comparison
does not mix the two-process iterative table with the one-process MOOSE run.

The repeated field-split regression uses two load steps on one reused PETSc
workspace. It prevents repeated solver setup from adding the same temperature
and mechanics index sets more than once. The full benchmark command is shown
below; replace `medium direct` with `medium hypre`, `large direct`, or
`large hypre` for the other measurements, and repeat each command once for the
warmed timing:

```bash
env \
  OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  MPIR_CVAR_CH4_NETMOD=ofi \
  FI_PROVIDER=tcp \
  /home/cooper/miniforge/envs/moose/bin/mpiexec -n 2 \
  taskset -c 0,1 \
  ./build/fuelsim_m1_single_core_benchmark medium direct 20 unscaled
```

These are measurements for two particular structured meshes, two MPI
processes, the current PETSc build, and the current default HYPRE subtype and
options. They are not a general strong-scaling result, a memory comparison, or
evidence that HYPRE will converge for another material, contact state, mesh, or
load path.

## 2026-08-11 M5.7 integrated transient parallel measurement

The M5.7 150-node, 104-element integrated transient case was measured at commit
`99d6e11` on a 13th-generation Intel Core i9-13980HX. PETSc 3.25.2 used direct
MUMPS factorization at every rank count. MPI ranks were restricted to CPUs
`0`, `0-1`, or `0-3`, respectively, and bound to cores. OpenMP, OpenBLAS, MKL,
and NumExpr were fixed to one thread. The benchmark executable writes no Exodus
or CSV result history; process wall time includes MPI launch, the solve, and the
small final flattened-state write or comparison.

The first run was kept separate. The warmed result is the median of three
successful runs with 97 accepted adaptive steps. Four ranks also produced one
valid 98-step run; it is reported separately instead of being mixed into the
equal-work statistic.

| MPI ranks | first wall time | warmed 97-step samples | warmed median | speedup over 1 rank | parallel efficiency |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 17.54 s | 17.16, 17.24, 17.23 s | 17.23 s | 1.000 | 100.0% |
| 2 | 19.99 s | 20.04, 20.05, 19.81 s | 20.04 s | 0.860 | 43.0% |
| 4 | 18.28 s | 18.26, 18.66, 18.32 s | 18.32 s | 0.941 | 23.5% |

Thus two ranks were `16.31%` slower and four ranks were `6.33%` slower than one
rank on the equal-step warmed statistic. The observed four-rank adaptive samples
`20.21`, `18.26`, and `18.66 s` have a median of `18.66 s`; the first sample used
98 accepted steps and the other two used 97. This small case therefore provides
no parallel speedup. Distributed assembly and MUMPS communication dominate the
work saved by partitioning 104 elements.

The same runs compared physical time, load factor, all three nodal fields, all
four components of elastic, plastic, and creep strain at every integration
point, both equivalent inelastic strains, all four stress components, contact
elastic slip, contact normal multiplier, and the stick/slip state. The largest
observed differences from the one-rank state were:

| path | temperature | displacement | strain history | stress | stick/slip state |
| --- | ---: | ---: | ---: | ---: | --- |
| 2 ranks, 97 steps | 5.3433e-11 K | 8.6225e-14 m | 6.1038e-11 | 1.4412 Pa | exact |
| 4 ranks, 97 steps | 5.3661e-11 K | 8.6052e-14 m | 6.1041e-11 | 1.4376 Pa | exact |
| 4 ranks, 98 steps | 3.2452e-4 K | 1.4628e-11 m | 2.9409e-9 | 242.26 Pa | exact |

All continuous contact histories also passed their absolute-plus-relative gates.
The 98-step differences are time-path differences caused by the iteration-based
adaptive threshold, not a failed distributed solve. They remain below the M5.7
limits of `1e-3 K`, `1e-10 m`, `1e-8` strain, and `1e3 Pa`. The one-rank state is
independently compared with the tracked MOOSE node, contact-pressure, material-
average, and 224-point integration data, so these complete-state comparisons
connect every measured rank count to the same MOOSE reference.

The one-rank reference command was:

```bash
env \
  OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  MPIR_CVAR_CH4_NETMOD=ofi \
  FI_PROVIDER=tcp \
  /usr/bin/time -f 'process_wall_seconds=%e' \
  taskset -c 0 \
  /home/cooper/miniforge/envs/moose/bin/mpiexec -bind-to core -n 1 \
  ./build/fuelsim_mpi_equivalence_benchmark \
  write_transient_integrated /tmp/m57_parallel_reference.txt \
  verification/fuelsim/transient_integrated_fuel_cladding.fsi \
  -pc_type lu -pc_factor_mat_solver_type mumps
```

For two ranks, the comparison command replaced the rank, CPU list, and mode:

```bash
env \
  OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  MPIR_CVAR_CH4_NETMOD=ofi \
  FI_PROVIDER=tcp \
  /usr/bin/time -f 'process_wall_seconds=%e' \
  taskset -c 0,1 \
  /home/cooper/miniforge/envs/moose/bin/mpiexec -bind-to core -n 2 \
  ./build/fuelsim_mpi_equivalence_benchmark \
  compare_transient_integrated /tmp/m57_parallel_reference.txt \
  verification/fuelsim/transient_integrated_fuel_cladding.fsi \
  -pc_type lu -pc_factor_mat_solver_type mumps
```

The four-rank command uses CPUs `0,1,2,3` and `-n 4`. These commands are manual
benchmarks and are intentionally not registered in the default CTest suite.
These measurements apply
only to this compact validation mesh, current adaptive path, hardware, PETSc,
and MUMPS build. They do not contradict the two-rank speedup measured on the
23,010-DOF benchmark and must not be generalized as a scaling limit.
