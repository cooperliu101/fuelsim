# Single-core benchmarks

`fuelsim_m1_single_core_benchmark` is the manual medium M1 benchmark. It uses:

```text
fuel mesh:       100 radial x 64 axial
cladding mesh:    16 radial x 64 axial
nodes/elements:  7,670 / 7,424
solution DOFs:   23,010
load steps:      20
linear solve:    PETSc LU (1 rank) or MUMPS (multiple ranks)
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
  ./build/fuelsim_m1_single_core_benchmark
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
