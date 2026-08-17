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

A follow-up lifted production problem dispatch out of the per-contribution
virtual call. Solver setup identifies the final steady or transient problem
once; callback loops use an explicit enum branch and direct calls, while custom
`NonlinearProblem` test fixtures retain the generic virtual path. With commit
`02309c7` as the immediate baseline, the same measurements were:

```text
                                      02309c7       candidate      observed change
default three-run median:            1.114839 s     1.108533 s       0.57% lower
medium 23,010 DOF, 20 steps:        60.722259 s    60.443432 s       0.46% lower
```

The small differences are consistent with eliminating one indirect call per
contribution, but the sample is only evidence of no regression rather than a
general speedup claim. Iteration counts and final residuals remained identical.

The next material-dispatch step encoded built-in thermal, elasticity, and
eigenstrain functions as enums when input functions were bound. Integration
points now use explicit switches and direct built-in calls; only user-registered
custom functions use the function-pointer fallback. With commit `562b00f` as
the immediate baseline, the default three-run median changed from `1.108533 s`
to `1.113615 s`, or 0.46% higher, while the single medium run changed from
`60.443432 s` to `60.408998 s`, or 0.06% lower. Iteration counts and residuals
were identical. These opposite sub-percent movements are treated as run-to-run
variation and no-regression evidence, not a speedup claim; the change is kept
for its explicit call graph, simpler debugging, and device-portability benefit.

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

## 2026-08-15 Cartesian three-dimensional inelasticity measurement

The Cartesian three-dimensional plasticity, creep, and coupled implementation was paired with its
pre-change commit `14ae880` in Release mode on CPU 0. OpenMP, OpenBLAS, MKL, and NumExpr were fixed to one
thread. The 1,584-DOF production input retained the first run only as a cold sample. Ten subsequent alternating
samples, including four pairs run in reverse order, had internal load-path medians of
`1.093807713/1.116627052 s` for the baseline/candidate. The candidate was `2.086%` slower. Both versions
completed 20 steps with 64 nonlinear and 64 linear iterations and the same `7.602509876370e-9` final residual.

The required 23,010-DOF, 20-step direct case completed once per version in `60.547771959/60.888084858 s`, a
`0.562%` candidate slowdown. Both runs used 62 nonlinear and 62 linear iterations, 82 residual callbacks, 62
Jacobian callbacks, one PETSc workspace, and the same `3.017657067419e-9` final residual. The matching MOOSE
run used the same mesh, physics, load steps, direct solver, CPU and thread limit, with CSV, Exodus and console
output disabled; its wall time was `47.32 s`. On this machine and current toolchain, MOOSE therefore used
`22.284%` less wall time than the candidate for this particular engineering-scale case. Internal fuelsim time
and MOOSE process wall time are not a general scaling comparison, but they establish that the new capability did
not change nonlinear work and that the remaining paired fuelsim slowdown is small rather than a large algorithmic
regression.

## 2026-08-14 source consolidation measurement

The source-consolidation candidate was paired with pre-refactor commit
`03b751e` in Release mode on CPU 0. OpenMP, OpenBLAS, MKL, and NumExpr were
fixed to one thread. The default 1,584-DOF direct case kept the first run
separate. Its following two internal load-path times had medians of
`1.116141079/1.163009159 s` for the baseline/candidate, so the candidate was
`4.20%` slower. Both versions completed 20 steps with 64 nonlinear and 64
linear iterations and the same `7.602509876370e-9` last residual.

The required 23,010-DOF, 20-step direct case was also completed once by each
version. The paired internal times were `60.769952996/61.888409483 s`, or a
`1.84%` candidate slowdown. Both used 62 nonlinear and 62 linear iterations,
82 residual callbacks, 62 Jacobian callbacks, one PETSc workspace, and the
same `3.017657067419e-9` last residual. These measurements establish numerical
work equivalence and a real performance regression in the original candidate.

The follow-up diagnosis reproduced a `7.22%` warmed slowdown over six paired
runs and localized the extra work to the consolidated steady/transient thermal
residual. The steady path passed zero-valued ADlite heat capacity and
temperature-rate scalars through the transient expression, so every thermal
node performed derivative-array multiplications that were absent before the
refactor. The shared function now uses null heat-capacity inputs to omit the
transient term without restoring duplicate steady and transient functions.

After that correction, seven alternating default-case runs kept the first pair
separate. The following six warmed samples had medians of
`1.091153993/1.097245903 s`; the corrected candidate was `0.56%` slower, within
the observed run-to-run spread. Two alternating engineering-case pairs gave
baseline/candidate means of `61.566550101/60.492010928 s`, so the corrected
candidate was `1.75%` faster. Mean Jacobian callback time changed from
`41.020744907` to `40.427390936 s`, or `1.45%` lower, and mean residual callback
time changed from `3.361947782` to `3.167908604 s`, or `5.77%` lower. Every run
retained the same iterations, callback counts, workspace count, and final
residual. These paired cases establish that the observed consolidation
regression was removed; they are not a broader performance or scaling claim.

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

## 2026-08-15 B3.3 three-dimensional contact regression check

The three-dimensional contact change was checked against pre-change commit
`ca14ec5`. Both Release builds used PETSc 3.25.2, one MPI process, CPU 0, and
one thread each for OpenMP, OpenBLAS, MKL, and NumExpr. One separate first run
was excluded before three interleaved warmed samples of the unchanged
1,584-degree-of-freedom RZ case:

| revision | warmed samples | median |
| --- | --- | ---: |
| `ca14ec5` | 1.082358, 1.081505, 1.087497 s | 1.082358 s |
| B3.3 candidate | 1.089770, 1.082828, 1.100855 s | 1.089770 s |

The candidate median is 0.68 percent higher. Both revisions completed exactly
64 nonlinear iterations, retained the same `7.602509876370e-09` final residual,
and created one PETSc workspace. This sub-percent movement is recorded as
no-regression evidence, not as a speedup or as a general performance claim.

The required 23,010-degree-of-freedom, 20-step fuelsim case then completed once
in `60.842814 s`, with 62 nonlinear iterations, 82 residual callbacks, 62
Jacobian callbacks, and one PETSc workspace. A warmed MOOSE run on the same
mesh and physics, with direct solve and file output disabled, took `46.35 s`.
Thus this fuelsim run was 31.3 percent slower than that MOOSE observation. The
comparison is a single engineering-case measurement and is not generalized to
other meshes or contact paths.

## 2026-08-15 dynamic contact assembly repair

The B3.3 production path kept every thermal and mechanical primary candidate as a runtime contribution. On the
23,010-degree-of-freedom RZ case this made the callback loop traverse mostly inactive zero blocks, assigned all volume
elements to rank 0 before contact candidates reached rank 1, and expanded the largest shadow state to all 23,010
degrees of freedom. The repair keeps the complete current-configuration primary-chain search and sparse-matrix graph,
but assembles only the selected thermal integration-point and mechanical-node candidates. PETSc's matrix preallocator
builds the complete graph once, and RZ candidate selection uses a double-only projection before the active ADlite
residual and Jacobian evaluation.

Release builds of pre-repair commit `a58c709` and the candidate used CPU 0, one MPI process, and one thread each for
OpenMP, OpenBLAS, MKL, and NumExpr. Separate first runs were excluded. The following three paired 1,584-degree-of-
freedom samples had internal load-path medians of `1.027653101/0.969833623 s`, so the candidate was 5.63 percent faster.
Both versions completed 64 nonlinear and linear iterations with the same `7.602509876370e-09` final residual.

The required 23,010-degree-of-freedom, 20-step direct case completed once per configuration:

| MPI ranks | `a58c709` | candidate | change |
| ---: | ---: | ---: | ---: |
| 1 | 60.842814 s | 30.528152 s | 49.82 percent lower |
| 2 | 28.682528 s | 12.084951 s | 57.87 percent lower |

Every medium run completed 62 nonlinear and linear iterations with 82 residual callbacks, 62 Jacobian callbacks, and
one PETSc workspace. On two ranks the largest shadow state fell from 23,010 to 12,197 degrees of freedom, and remote
shadow values fell from 11,672 to 8,057. The previously recorded matching output-disabled MOOSE warmed time is
`46.35 s`; the repaired one-rank fuelsim observation is 34.14 percent lower. The repaired two-rank path is still 17.8
percent above the older `a462afe` specialized-path observation of `10.260384 s`, because production contact now
searches the complete current primary chain and supports large sliding instead of assuming a fixed near-neighbor
projection. These measurements apply only to the named mesh, physics, direct solver, and CPU placement.

## 2026-08-15 exact large-surface contact search

The current-configuration contact search now caches the previous active primary segment or face and uses a refitted
axis-aligned bounding-box tree when a contact surface has more than 64 primary candidates. The tree only prunes a node
when its Euclidean distance lower bound is strictly greater than the best exact projected distance, so equal-distance
candidates still reach the existing deterministic primary-index tie break. Smaller surfaces retain the compact linear
scan because paired measurements showed that building and traversing a tree for the 64-segment engineering case cost
more than it saved. The sparse matrix still reserves every potential candidate.

A deterministic 97-segment unit check compared every tree result with exhaustive selection over 151 query points. It
visited at most one exact candidate per query while returning the same nearest item, and a separate equal-distance case
retained both candidates for the primary-index tie break. The 128-segment M5.2 large-sliding solve retained its
cross-segment ownership and MOOSE field errors; direct process observations were `0.5080/0.5073 s` for the candidate and
`dabfcb2`, which is indistinguishable at this scale rather than evidence of an end-to-end speedup.

All final timings used CPU 0 and one thread for OpenMP, OpenBLAS, MKL, and NumExpr. Three final 1,584-degree-of-freedom
candidate samples were `1.041394112`, `1.039736670`, and `1.067377675 s`, with a `1.041394112 s` median. Three
`dabfcb2` samples were `1.039641740`, `1.067660264`, and `1.040072855 s`, with a `1.040072855 s` median. The candidate
median is 0.13 percent higher, which is recorded as no material change rather than a speedup. The required final
23,010-degree-of-freedom, 20-step run completed in `31.538306089 s`; the two same-session `dabfcb2` observations were
`31.942540778` and `31.708051782 s`. The final two-rank candidate and paired baseline observations were
`11.989104812/11.923955247 s`. Every medium run completed 62 nonlinear and linear iterations, 82 residual callbacks, 62
Jacobian callbacks, and one PETSc workspace. These results show that the search changes preserve current benchmark
performance; the present 64-segment benchmark does not exercise the tree and therefore does not establish large-surface
end-to-end speedup.

## 2026-08-16 M5.8 cylindrical PCMI Hex8 comparison

M5.8 is now a uniformly refined three-dimensional benchmark. Its solid fuel
uses 640 Hex8 elements and its independent annular cladding uses 512, for a
total of 1,617 nodes, 1,152 elements, and 6,468 coupled degrees of freedom. The
fuel surface has 144 mechanical contact nodes and 128 Quad4 faces; the cladding
inner surface has 256 Quad4 faces. The case advances 20 fixed Backward Euler
steps. Both programs read the same tracked Exodus mesh and use the same
finite-strain material laws, heat and pressure histories, penalty contact,
friction coefficient, direct LU solve, and nonlinear tolerances.

CPU 0 was fixed with `taskset`, while OMP, OpenBLAS, MKL, and NumExpr were each
limited to one thread. Field and console output were disabled for the recorded
timings. One complete fuelsim run and the MOOSE reference-generation run were
used as unrecorded warm-ups. The subsequent fuelsim wall clock was `343.47 s`;
the subsequent MOOSE wall clock was `777.90 s`. The measured ratio is `2.2648`,
or a 55.85 percent lower fuelsim wall time. The fuelsim comparison run used 75
nonlinear iterations and one PETSc workspace. Because each sample costs several
minutes, these are single warmed observations rather than three-run medians.
They establish performance for this exact one-process direct-solver benchmark,
not parallel scaling or engineering-size asymptotic behavior.

From the repository root, the recorded fuelsim command was:

```bash
env PATH=/home/cooper/miniforge/envs/moose/bin:/usr/local/bin:/usr/bin:/bin \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 taskset -c 0 \
  /usr/bin/time -f 'wall_seconds=%e' \
  ./build/fuelsim -i verification/fuelsim/transient_integrated_hex8.fsi
```

From `verification/moose`, the matching MOOSE command was:

```bash
env PATH=/home/cooper/miniforge/envs/moose/bin:/usr/local/bin:/usr/bin:/bin \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 taskset -c 0 \
  /usr/bin/time -f 'wall_seconds=%e' \
  /home/cooper/projects/july/july-opt -i m58_integrated_hex8.i \
  Outputs/csv=false Outputs/console=false
```

The final MOOSE field comparison uses all 1,617 nodes, all 144 contact nodes,
and all 512 cladding elements. Temperature, three displacement components,
contact pressure, equivalent stress, equivalent plastic strain, and equivalent
creep strain each pass relative L2, relative absolute-peak, and maximum
pointwise relative errors below 0.5 percent without a denominator floor. The
largest default acceptance metric is the X-displacement maximum pointwise
relative error at `0.010724 percent`.

The full field comparison is intentionally not registered with CTest. Run it
manually from the repository root with:

```bash
./build/fuelsim_m58_integrated_hex8_benchmark \
  verification/fuelsim/transient_integrated_hex8.fsi \
  verification/moose/m58_integrated_hex8_all_nodes_0020.csv \
  verification/moose/m58_integrated_hex8_contact_pressure_0020.csv \
  verification/moose/m58_integrated_hex8_clad_state_0020.csv
```

## 2026-08-17 latest-code 6,468-DOF one-, two-, and four-process comparison

The latest fuelsim executable at commit `623a67a` was paired with the same
tracked M5.8 Exodus mesh and direct MUMPS input at one, two, and four MPI
processes. Every fuelsim and MOOSE run completed all 20 fixed time steps. CSV
and Exodus output were disabled; MOOSE console output was enabled only so that
long runs could be monitored with line-buffered progress output.

| MPI processes | fuelsim wall time | MOOSE wall time | fuelsim speedup over MOOSE | fuelsim parallel efficiency | MOOSE parallel efficiency |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `78.23 s` | `777.19 s` | `9.9347x` | `100.00%` | `100.00%` |
| 2 | `45.58 s` | `489.10 s` | `10.7306x` | `85.82%` | `79.45%` |
| 4 | `30.44 s` | `270.72 s` | `8.8936x` | `64.25%` | `71.77%` |

Relative to the corresponding MOOSE process count, fuelsim was 89.93 percent,
90.68 percent, and 88.76 percent faster at one, two, and four processes. The
fuelsim four-process result is just below the 65 percent efficiency target
under this latest single-sample wall-clock measurement; it should not be
confused with the separate 30,148-DOF GMRES field-split result. These are
single warmed observations on the named mesh, hardware, PETSc/MUMPS build, and
MOOSE executable, not general scaling limits.

## 2026-08-16 M5.8 two-process parallel efficiency

The original runtime contribution partition divided the 2,688 active blocks by
count. Its two intervals were `[0,1344)` and `[1344,2688)`. Because the first
1,152 blocks are eight-point finite-strain Hex8 volume integrations, rank zero
performed all volume work while rank one mostly handled the cheaper surface
blocks. With one-process and two-process MUMPS runs using the same 76 nonlinear
iterations, the measured internal times were `111.122 s` and `103.712 s`; the
two-process efficiency was only `53.57 percent`.

The runtime partition now balances operation-based work estimates while retaining
one contiguous and unique interval per process. Thermal contact, mechanical
contact, and pressure Jacobians submit only their current nonzero row and column
subblocks to PETSc; the complete preallocated graph is unchanged. The Cartesian
three-dimensional sparsity graph also stores one representative for each
secondary-face/primary-face pair instead of repeating the same 32-DOF block for
four thermal quadrature points and four mechanical face nodes. Runtime physics,
candidate search, and every active contribution remain unchanged.

The final two-process MOOSE comparison used intervals `[0,664)` and `[664,2688)`.
They are contiguous, nonoverlapping, and cover every active contribution. The
largest shadow state was 3,972 DOFs; both shadows contained 7,751 values in total,
including 3,672 remote values. The complete solve used 75 nonlinear iterations,
`1.716526 s` of setup, `6.020732 s` in residual callbacks,
`35.818805 s` in Jacobian callbacks, and `59.434407 s` internally. All M5.8
nodal, contact-pressure, and material-state comparisons remained below the
`0.5 percent` three-metric limit.

The final end-to-end production-entry measurements used a newly linked Release
executable, the same MUMPS direct solver for both process counts, CPUs 0 and 1,
and one thread for every numerical library. Three samples were collected for
each process count; their wall times were `104.37`, `104.03`, and `103.83 s` for
one process, and `64.43`, `64.87`, and `65.40 s` for two processes:

| MPI processes | median process wall time | speedup | parallel efficiency |
| ---: | ---: | ---: | ---: |
| 1 | `104.03 s` | `1.0000` | `100.00 percent` |
| 2 | `64.87 s` | `1.6037` | `80.1834 percent` |

Thus the median two-process wall time is `37.64 percent` lower and its parallel
efficiency exceeds the `80 percent` target for this exact 6,468-DOF, 20-step
benchmark. This is a two-process strong-scaling result for the named mesh,
physics, direct solver, hardware placement, and toolchain; it is not evidence of
general multi-node or larger-process-count scaling.

The one-process command was:

```bash
env PATH=/home/cooper/miniforge/envs/moose/bin:/usr/local/bin:/usr/bin:/bin \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 MPIR_CVAR_CH4_NETMOD=ofi FI_PROVIDER=tcp \
  PETSC_OPTIONS='-pc_factor_mat_solver_type mumps' \
  /usr/bin/time -f 'process_wall_seconds=%e' \
  taskset -c 0 \
  /home/cooper/miniforge/envs/moose/bin/mpiexec -bind-to core -n 1 \
  ./build/fuelsim -i verification/fuelsim/transient_integrated_hex8.fsi
```

The tracked M5.8 input now explicitly selects MUMPS at every process count, so
`PETSC_OPTIONS` is no longer required. The two-process command changes the CPU
list and process count to `taskset -c 0,1` and `-n 2`.

## 2026-08-16 M5.8 four-process direct-MUMPS efficiency

The four-process work keeps the direct MUMPS algorithm. The M5.8 input selects
MUMPS for one process as well as multiple processes, uses the same SCOTCH
ordering unless the user overrides it, and retains the four-process MUMPS
memory-relaxation guard. PETSc uses an internal node-major ordering to keep the
four fields of each Hex8 node local while the public state remains field-major.
Exact local Jacobian patterns omit structurally zero thermal-mechanical blocks,
and the four contiguous contribution intervals are `[0,346)`, `[346,692)`,
`[692,1038)`, and `[1038,2688)`.

Three complete one-process and four-process MOOSE-comparison runs used the same
20 fixed time steps, 93 nonlinear iterations, 49 Jacobian evaluations, direct
MUMPS factorization, and one thread per numerical library. The internal total
solver times were `73.479816`, `73.504522`, and `73.741511 s` for one process,
and `24.965830`, `25.082956`, and `24.926659 s` for four processes. Their
medians give speedup `2.944205` and four-core solver efficiency `73.6051
percent`, above the 70 percent target. The nonlinear-solve-only median efficiency
is `75.6576 percent`.

The corresponding process wall times were `77.85`, `77.86`, and `78.12 s`, and
`30.29`, `30.41`, and `30.21 s`. Their medians give `64.2621 percent`
end-to-end efficiency because replicated Exodus input and problem construction
remain outside the solver timer. This fixed overhead is reported separately and
is not claimed to exceed 70 percent.

Every run passed the full MOOSE comparison. A one-process reference followed by
four-process degree-of-freedom comparison found at most `4.43e-10 K`
temperature difference and `2.52e-12 m` displacement difference; the largest
configured tolerance ratio was `0.02491`. The four intervals are contiguous,
nonoverlapping, and cover all 2,688 contributions. These measurements establish
strong scaling only for this named 6,468-DOF benchmark and hardware.

## 2026-08-16 compact three-dimensional contact candidates

The 30,148-degree-of-freedom M5.8 extension has 640 secondary fuel faces and
1,280 primary cladding faces. The previous Cartesian contact construction stored
3,276,800 complete thermal candidates and 3,276,800 complete mechanical
candidates, including repeated coordinates, shape functions, derivatives, and
node arrays. The compact implementation stores each primary and secondary face
once, represents the complete candidate product by offsets and indices, and
constructs only candidates visited by the exact search or active assembly. The
PETSc sparsity preallocation still visits every secondary-face/primary-face pair,
so the admissible contact graph, unique projection, and large-sliding behavior do
not change.

One process on CPU 0 with every numerical library fixed to one thread completed
the same first `0.05 s` step before and after the change. Both runs used 30,148
DOFs, six nonlinear and linear iterations, and one PETSc workspace. After
excluding memory and timing fields, their complete console outputs were
identical, including residuals, conservation diagnostics, material summaries,
and contact results.

| measurement | `77b1dcc` | compact candidates | change |
| --- | ---: | ---: | ---: |
| initial resident memory | 3.45 GiB | 114.43 MiB | 96.76 percent lower |
| PETSc maximum resident memory | 5.41 GiB | 1.73 GiB | 67.93 percent lower |
| Linux process high-water mark | 5.84 GiB | 2.01 GiB | 65.66 percent lower |
| internal total time | 92.1779 s | 82.7908 s | 10.18 percent lower |
| monitored process wall time | 100 s | 86 s | 14.00 percent lower |

The unrelated 1,584-DOF RZ path was paired against `77b1dcc` after one warm-up.
Three interleaved internal-time samples had medians of `1.016027090 s` and
`1.025738838 s`, a 0.96 percent compact-candidate time increase with identical 63 nonlinear
and linear iterations and the same `8.362012981744e-09` final residual. The
required 23,010-DOF, 20-step RZ case also completed in `28.911380536 s` with 62
nonlinear and linear iterations, 82 residual callbacks, 62 Jacobian callbacks,
and one PETSc workspace. The 30,148-DOF result covers one time step and is an
exact memory and paired-speed observation for this mesh, not a general scaling
claim.

## 2026-08-17 M5.8 30,148-DOF four-process efficiency

The original sparse preallocation divided contact candidates by candidate
index. Each process therefore inserted most of its preallocation rows into
other processes' matrix ownership ranges. On four processes, PETSc spent
`60.112 s` in the beginning of matrix assembly before the nonlinear solve. The
one-process and four-process direct-MUMPS internal times were `82.901277 s` and
`104.420422 s`, respectively, so the original four-core efficiency was only
`19.85 percent`.

The corrected preallocation visits the complete deterministic graph on every
process but inserts only locally owned rows. It does not change runtime
contribution ownership, residual values, Jacobian values, or the matrix graph.
This reduced the four-process preallocation assembly begin from `60.112 s` to
`0.641 s`; the otherwise unchanged four-process direct-MUMPS time fell to
`42.541975 s`. SCOTCH remained the fastest tested MUMPS ordering: the same
one-step direct case took `42.42 s` with SCOTCH, `47.27 s` with PT-SCOTCH,
`47.43 s` with METIS, and `68.48 s` with PORD.

The final 30,148-DOF input uses GMRES with multiplicative temperature and
mechanics field splitting. Here GMRES is the Krylov linear solver, and the field
split is a physics-based preconditioner. On four or more processes, the
mechanics block defaults to one level of incomplete-LU fill unless the user
explicitly provides another PETSc option. Zero fill was faster for the first
step but was rejected: the complete 20-step path stopped after 14 accepted
steps with a divergent linear solve. One fill level completed all 20 steps.

One fixed `0.05 s` step used CPUs 0 through 3, the default shared-memory MPI
transport, and one thread for OpenMP, OpenBLAS, MKL, and NumExpr. No TCP
fallback transport was forced. The paired results were:

| measurement | one process | four processes | speedup | four-core efficiency |
| --- | ---: | ---: | ---: | ---: |
| internal total time | 81.708530 s | 21.385841 s | 3.8207 | 95.52 percent |
| process wall time | 82.51 s | 23.18 s | 3.5595 | 88.99 percent |
| linear iterations | 482 | 732 | - | - |

The four-process solution was also compared degree by degree with the
one-process direct-MUMPS reference. The maximum temperature difference was
`1.59844e-10 K`, the maximum displacement difference was `3.84966e-13 m`, and
the largest configured tolerance ratio was `0.0038491`. The four contribution
intervals were contiguous, nonoverlapping, and covered all contributions.

The selected one-fill configuration completed the full 20-step four-process
path before it was installed as the default: 20 steps were accepted with no
cutback or rejected step, using 138 nonlinear iterations, 16,882 linear
iterations, one PETSc workspace, `325.087564 s` internal time, and `326.91 s`
wall time. The aggregate process high-water mark was `1,388,904,448 bytes`
(`1.29 GiB`), and the largest individual process high-water mark was
`356,122,624 bytes`. A subsequent one-step run without an explicit fill option
reproduced 732 linear iterations and the timing above, confirming that the
program default selects the tested configuration.

The unrelated 1,584-DOF RZ path was paired with pre-change commit `6b7c995` on
CPU 0. After one warm-up, the three-run internal-time medians were
`0.991434761 s` before and `0.991832096 s` after, a `0.04 percent` increase;
both completed 63 nonlinear iterations. The required current 23,010-DOF,
20-step direct case completed in `28.445658 s` with 62 nonlinear and linear
iterations and one PETSc workspace. The matching output-disabled MOOSE run
took `46.77 s` wall time. Finally, all 60 CTest registrations passed serially.
These results establish greater than 65 percent strong-scaling efficiency only
for the named 30,148-DOF case, first-step workload, CPU placement, PETSc build,
and solver configuration.
