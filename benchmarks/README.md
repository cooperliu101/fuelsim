# Solver benchmarks

2026-09-09: Axisymmetric MOOSE comparison inputs and the dedicated PCMI parallel
comparison script were removed. Commands and measurements below that refer to
those files are historical records; they are not current runnable comparisons.

2026-09-09：原来的 C++ 轴对称性能入口已删除。当前使用完整的生产输入卡：

```bash
build/fuelsim -i verification/fuelsim/steady_rz_performance_medium.fsi
build/fuelsim -i verification/fuelsim/steady_rz_performance_large.fsi
```

两档分别为 23,010 和 45,630 个自由度，采用 CAX4T、小应变、稳态热弹性接触，
固定 20 个加载增量。与原程序的默认 Quad4 单元公式不同，因此下面旧记录的
速度不能直接用于新输入。精度与速度的当前证据见
[轴对称生产性能验证](../verification/abaqus/rz_performance/README.md)。
同名 `_timing.fsi` 输入卡关闭结果文件输出，用于独立计时。
`run_rz_performance.py` 顺序运行两套程序，每套先预热一次，再记录两次；
它不会生成或修改输入卡。`--size medium` 可以只运行中等规模，
`--results-directory` 和 `--windows-results` 指定新的证据目录。中等规模完整牛顿更新
结合载荷预测，在 0.01% 精度要求下，非线性迭代从 65 次降到 46 次；
ADlite 0.2.3 下外部时间均值为 22.44 秒，Abaqus 同批次为 41.82 秒；
本次旧版 ADlite 0.2.2 重测为 22.83 秒，更新后的改善约为 1.68%；
大规模仍保留原配置。下文原程序的命令、预条件器研究和性能记录均属于历史资料。

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
cmake --build build --parallel 4

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
  ./build/fuelsim -i verification/fuelsim/transient_integrated_c3d8t.fsi
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
  verification/fuelsim/transient_integrated_c3d8t.fsi \
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

The latest one-process full-field comparison also passed using the tracked MOOSE
reference files. It completed 20 steps with 93 nonlinear iterations, 49
Jacobian evaluations, and one PETSc workspace in `75.9478 s` internal time
(`80.38 s` process wall time). The largest maximum pointwise relative error was
`1.07366e-04` (`0.0107366 percent`) for X displacement; temperature,
displacements, contact pressure, equivalent stress, equivalent plastic strain,
and equivalent creep strain all remained below the `0.5 percent` three-metric
limit. The final state had 144 active contact nodes, 109 sliding contact nodes,
and 144 nodes crossing primary faces.

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
  ./build/fuelsim -i verification/fuelsim/transient_integrated_c3d8t.fsi
```

The tracked M5.8 input now explicitly selects MUMPS at every process count, so
`PETSC_OPTIONS` is no longer required. The two-process command changes the CPU
list and process count to `taskset -c 0,1` and `-n 2`.

## 2026-08-16 M5.8 four-process direct-MUMPS efficiency

This historical four-process comparison used the direct MUMPS algorithm. The
M5.8 input selected MUMPS for one process as well as multiple processes and used
the then-default SCOTCH ordering. The code retained the four-process MUMPS
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

The final 30,148-DOF input now defaults to direct MUMPS. The measurements below
used the earlier GMRES configuration with multiplicative temperature and
mechanics field splitting; GMRES is the Krylov linear solver, and the field
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

## 2026-08-19 Hex8 narrow two-level constitutive AD

The three-dimensional Hex8 Jacobian path no longer seeds all 32 local DOFs with
ADlite. The kinematics chain is seeded on the nine displacement-gradient
components plus the quadrature-point temperature (width 10); the constitutive
evaluation is seeded on the six strain components plus temperature (width 7) and
reattached to the kinematics chain with `adlite::compose`; the final chain from
the point seeds to the 32 local DOFs is linear and applied in closed form. The
residual-only path, the RZ kernels, the Quad4Face boundary, and the
three-dimensional contact kernels are unchanged, and the assembled Jacobian
remains the exact consistent tangent: the local centered directional-difference
checks and the complete CTest suite pass unchanged.

All runs used Release builds, one MPI process pinned to CPU 0, and one thread
for OpenMP, OpenBLAS, MKL, and NumExpr. The paired baseline is the pre-change
commit `ad4cdf8`. One warm-up run preceded the recorded samples on each side.

M5.8 6,468-DOF, 20-step finite-strain Hex8 case
(`fuelsim_m58_integrated_hex8_benchmark` with the three tracked MOOSE CSV
references), two recorded samples per side:

| measurement | `ad4cdf8` baseline | narrow-AD candidate | change |
| --- | ---: | ---: | ---: |
| Jacobian callback seconds | 46.871 / 48.291 | 25.452 / 25.447 | about 46.5 percent lower |
| residual callback seconds | 12.211 / 15.209 | 12.190 / 12.180 | unchanged |
| process wall seconds | 77.76 / 82.66 | 55.98 / 55.96 | about 29.4 percent lower |

Every run completed 20 steps with 93 nonlinear iterations, 113 residual
evaluations, 49 Jacobian evaluations, and one PETSc workspace. The full-field
MOOSE comparison kept every field at its pre-change error magnitude: the largest
maximum pointwise relative error stayed at `1.0737e-4` for X displacement, and
temperature, all displacement components, contact pressure, equivalent stress,
and both equivalent inelastic strains remained below the 0.5 percent
three-metric limit.

The untouched RZ path showed no regression. The default 1,584-DOF steady case
(`fuelsim -i verification/fuelsim/steady_fuel_cladding.fsi`) kept 20 load steps,
63 nonlinear iterations, one workspace, and the same `8.362012981744e-9` final
residual on both sides; warmed three-run medians were `0.950063 s` for the
baseline and `0.944727 s` for the candidate, a 0.56 percent difference within
the observed run-to-run spread. The 23,010-DOF, 20-step medium case completed
once per side in `28.917 s` (baseline) and `28.241 s` (candidate) with identical
62 nonlinear iterations, 82 residual and 62 Jacobian callbacks, one workspace,
and the same final residual `3.374857832964e-9`. These RZ differences are
run-to-run variation on this machine and are recorded as no-regression evidence,
not as a speedup claim.

## 2026-08-19 RZ narrow two-level constitutive AD

The two-dimensional axisymmetric RZ Quad4 Jacobian path now uses the same
narrow two-level seeding as the Hex8 path. The kinematics chain is seeded on
the four in-plane displacement-gradient components plus the quadrature-point
radial displacement and temperature (width 6); the constitutive evaluation is
seeded on the four strain components `[rr, zz, hoop, rz]` plus temperature
(width 5) and reattached to the kinematics chain with `adlite::compose`; the
final chain from the point seeds to the 12 local DOFs is linear and applied in
closed form. The residual-only path is unchanged, the assembled Jacobian
remains the exact consistent tangent, and the residual values are bitwise
identical to the pre-change build; only the floating-point summation order of
the Jacobian entries changed. The local centered directional-difference checks
and the complete CTest suite pass.

All runs used Release builds, one MPI process pinned to CPU 0, and one thread
for OpenMP, OpenBLAS, MKL, and NumExpr. The paired baseline is the pre-change
commit `9067649`. One warm-up run preceded the recorded samples on each side.

Default 1,584-DOF steady case (`fuelsim -i
verification/fuelsim/steady_fuel_cladding.fsi`), three recorded samples per
side, medians:

| measurement | `9067649` baseline | narrow-AD candidate | change |
| --- | ---: | ---: | ---: |
| internal total seconds | 1.029771 | 0.964532 | about 6.3 percent lower |
| nonlinear iterations | 63 | 64 | one additional iteration |
| final residual norm | 8.362012981744e-9 | 8.316025877667e-9 | last-digit level |

The small case takes one additional nonlinear iteration on the candidate. That
is the expected consequence of the changed Jacobian summation order under the
iteration-count-sensitive SNES path, not a physics change: with the case's
`target_nonlinear_iterations` step-size control disabled, the pre-change and
post-change builds produce bitwise-identical load paths.

23,010-DOF, 20-step medium case (`fuelsim_m1_single_core_benchmark medium
direct 20 unscaled`), one run per side:

| measurement | `9067649` baseline | narrow-AD candidate | change |
| --- | ---: | ---: | ---: |
| Jacobian callback seconds | 12.614745 | 10.340405 | about 18.0 percent lower |
| residual callback seconds | 2.748301 | 2.717152 | about 1.1 percent lower |
| nonlinear solve seconds | 28.877892 | 26.438876 | about 8.4 percent lower |
| load-path total seconds | 29.164330 | 26.715794 | about 8.4 percent lower |

Both sides completed 20 steps with 62 nonlinear and linear iterations, 82
residual and 62 Jacobian callbacks, and one PETSc workspace. The final
residual differed only at the level of the summation-order rearrangement
(`3.374857832964e-9` versus `3.249125466882e-9`).

The M5.7 integrated validation case is sensitive to the same iteration-count
effect: its step-size controller targets 8 plus-or-minus 2 nonlinear
iterations, so the candidate's accepted-step count changed from 26 to 18 while
the solution itself did not move (with the controller's iteration target
disabled, both builds produce the identical 17-step path and agree to twelve
significant digits in the maximum error estimate). The case's `time_sequence`
and the five tracked MOOSE CSV references were regenerated for the 18-step
path with the recorded `july-opt` binary; every M5.7 metric remains far below
its 0.5 percent limit (largest field metric 0.206719 percent, largest
quadrature-point metric 0.148170 percent).

## 2026-08-30 C3D8RT contact timing and two-process equivalence

The C3D8RT baseline reuses the M5.8 cylindrical 1,617-node, 1,152-element,
6,468-degree-of-freedom mesh and its 20 fixed Backward Euler steps. Both fuel
and cladding regions select C3D8RT in
`verification/fuelsim/transient_integrated_c3d8rt.fsi`; thermal and mechanical
contact, Coulomb friction, cladding plasticity and creep, and one PETSc
workspace remain active. OpenMP and library thread counts were one. These runs
used the same executable and environment but were not CPU-pinned, so the
numbers are a baseline and cost-location measurement rather than a general
speed claim.

The original C3D8RT implementation evaluated the complete finite-strain
residual once for each of 32 Jacobian columns. It now preserves the width-7
material and width-10 kinematics chains, then evaluates one 24-displacement
geometry block and one eight-temperature block. It does not use a complete
32-degree-of-freedom identity seed. The residual, material data, and fixed
`0.005` hourglass coefficient are unchanged.

The same CTest command before and after this change gives:

| measurement | column-by-column C3D8RT | two-block C3D8RT | change |
| --- | ---: | ---: | ---: |
| wall seconds | 452.58 | 135.10 | 70.15 percent lower |
| Jacobian callback seconds | 379.651 | 54.6117 | 85.62 percent lower |

A separate explicitly one-thread sequential pair used identical commands
apart from the element choice:

| measurement | C3D8T, one process | optimized C3D8RT, one process |
| --- | ---: | ---: |
| wall seconds | 203.38 | 127.59 |
| residual callback seconds | 84.6926 | 40.8294 |
| Jacobian callback seconds | 79.9235 | 53.3877 |

On this fixed case, optimized C3D8RT wall time is `0.62734` times C3D8T, or
`37.27%` lower, and its Jacobian callback time is `33.20%` lower. These are
paired measurements for this machine and path, not a general element-speed
claim.

The release-only CTest fixture first writes the complete one-process
degree-of-freedom vector and then compares the two-process result. The one-
process contribution interval is `[0,2688)`. The two-process intervals are
`[0,664)` and `[664,2688)`, which are nonoverlapping and complete. Maximum
one-process/two-process differences are `1.393800630467e-9 K` for temperature
and `6.565067039261e-12 m` for displacement; the largest difference divided by
its explicit tolerance is `0.06513751675647`. The optimized fixture takes
`135.10 s` at one process and `79.20 s` at two processes, or `214.31 s` total.
This is direct evidence for the
C3D8RT region-plus-contact combination rather than merely evidence that two
processes can start.

`verification/fuelsim/transient_integrated_c3d8rt_30k.fsi` also records the
30,148-degree-of-freedom manual C3D8RT configuration. It was not timed in this
change; the larger configuration is therefore a reproducible pending
benchmark, not a reported result.

### 2026-08-30 C3D8RT finite-Jacobian follow-up

The paired baseline is commit `53378b5`. Both Release runs used one MPI
process pinned to CPU 0 and one thread for OpenMP, OpenBLAS, MKL, and NumExpr.
The CPU governor reported `performance`, but hardware frequency was not fixed.
One unrecorded baseline run warmed the executable and input path before the
recorded pair. Both sides solved the same 6,468-degree-of-freedom, 20-step
`transient_integrated_c3d8rt.fsi` case and wrote the complete degree-of-freedom
vector.

The follow-up removes two remaining sources of repeated work. A Jacobian call
now takes its primal residual from the 24-displacement block instead of first
performing a separate passive residual evaluation. The eight-temperature block
reuses the already available passive current geometry and central displacement
gradient, and evaluates only terms with nonzero temperature derivatives. It
therefore omits only the temperature-independent body-source and mechanical
hourglass terms. The fixed `0.005` hourglass coefficient, quadrature, material
properties, residual, and complete consistent Jacobian are unchanged.

| measurement | `53378b5` baseline | follow-up candidate | change |
| --- | ---: | ---: | ---: |
| end-to-end wall seconds | 130.13 | 119.37 | 8.27 percent lower |
| residual callback seconds | 40.8384 | 40.9964 | 0.39 percent higher |
| Jacobian callback seconds | 55.4484 | 44.2830 | 20.14 percent lower |

The baseline and candidate degree-of-freedom files are both 132,341 bytes and
have the identical SHA-256 digest
`e56452c4e5c0cd3ca6c6367d1d373898e6de9d425e90c81c61f3a20afe6cb558`.
The local centered-difference tests also require the residual returned by a
Jacobian call to equal the residual-only path entry by entry.

## M5.8 C3D8T Abaqus timing

The B5.46 Abaqus comparison converts the exact M5.8 Exodus mesh with
`verification/abaqus/exodus_to_abaqus.py`. Abaqus R2018x uses one process,
twenty fixed increments, and final-only field output. Its external wall time is
`36.115237 s`; the analysis summary reports `24.8 s` CPU and `26 s` wall time.
The final pinned, one-thread Fuelsim C3D8T Release external-wall trials take
`34.15 s`, `34.08 s`, and `33.99 s` when both Fuelsim and ADlite use
interprocedural optimization. Their `34.08 s` median is `5.635%` below Abaqus
and gives an observed ratio of `0.943646`. The representative full-field Abaqus
comparison takes `34.35 s` externally and reports `27.626988 s` internally, of
which residual and Jacobian callbacks consume `6.487015 s` and `12.795092 s`.
It uses 92 nonlinear iterations, 112 residual evaluations, and 26 Jacobian
evaluations, versus Abaqus's 72 iterations and decompositions. The C3D8T
residual-only path reuses ordinary-double finite-kinematics and stress caches
instead of constructing the derivative-bearing element-pressure system, while
matching the active path's arithmetic order. A Jacobian reuse period of five
reduces derivative assembly without changing the fixed time steps or convergence
limits. The comparison crosses native Windows and WSL and is not treated as a
pure kernel comparison. No material, contact, load, convergence, or
stabilization coefficient changes.

The required unrelated RZ non-regression paths used the same pinning and
thread settings. After one warm-up per executable, the 1,584-degree-of-freedom
pre-change internal times were `0.861824 s`, `0.858642 s`, and `0.854164 s`;
the candidate times were `0.852278 s`, `0.849842 s`, and `0.854755 s`. The
medians are `0.858642 s` and `0.852278 s`, respectively, or 0.74 percent lower
for the candidate. The 23,010-degree-of-freedom, 20-step candidate completed in
`25.355092 s` with 62 nonlinear iterations, 82 residual evaluations, 62
Jacobian evaluations, one PETSc workspace, and a final residual norm of
`3.249125466882e-09`.

## M5.8 C3D8RT Abaqus timing

B5.47 converts the exact same 1,617-node and 1,152-element M5.8 Exodus mesh to
Abaqus C3D8RT. Both Abaqus and Fuelsim complete twenty fixed `0.05 s`
increments without cutback. The controlled one-process Abaqus R2018x run takes
`23.489094 s` external wall time and 74 iterations. The paired pre-change
Fuelsim trial at commit `6cf276c` takes `129.70 s`. After the passive residual,
closed geometry/Jacobian chain, explicit backtracking line search, and a
Jacobian reuse period of three, the predictor-off CPU-0-pinned Fuelsim Release
control takes `29.99 s`. Enabling the recent-two-committed-step linear predictor
reduces three final repeated trials to `22.32 s`, `22.44 s`, and `22.75 s`; their
`22.44 s` median is `25.18%` below the paired predictor-off control and `4.47%`
below Abaqus. The representative `18.062041 s` internal total contains
`2.381772 s` in residual callbacks and `4.611456 s` in Jacobian callbacks, with
84 nonlinear iterations, 104 residual evaluations, and 35 Jacobian evaluations.
The Fuelsim-to-Abaqus external-wall ratio is `0.955336975`. This is a controlled
end-to-end engineering comparison across native Windows and WSL, not a general
speed ratio or a pure element-kernel or linear-solver comparison.

The fixed-path field qualification is documented by B5.47 rather than inferred
from timing. All ordinary scalar and axial field metrics pass `0.5%`; the
radial pointwise error uses the recorded `4%` small-reference gate, complete
Cartesian and tangential diagnostics remain visible, and analytical-zero
tangential displacement uses a `1 um` absolute bound. Abaqus artificial strain
energy stays below `0.221711%` of internal energy. No physical or hourglass
coefficient is changed for either the error or timing comparison. The predictor
changes only the Newton initial guess; its previous committed node state and
time are included in rollback, state snapshots, and checkpoint version 17.

## M5.8 finite-sliding surface-to-surface contact timing

B5.53 and B5.54 repeat the same 1,617-node, 1,152-element, 6,468-degree-of-
freedom, twenty-increment paths with explicit finite-sliding surface-to-surface
mechanical contact in both Fuelsim and Abaqus. The element types are C3D8T and
C3D8RT, respectively. Materials, loads, friction, penalty, convergence limits,
time increments, direct MUMPS solve, and final-only reference output are
unchanged from B5.46 and B5.47.

CPU 0 is pinned and `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
`MKL_NUM_THREADS`, and `NUMEXPR_NUM_THREADS` are all one. The controlled
single-sample external-wall results are:

| Case | Fuelsim | Abaqus | Fuelsim reduction | Fuelsim/Abaqus |
| --- | ---: | ---: | ---: | ---: |
| C3D8T B5.53 | `33.52 s` | `35.480477 s` | `5.5255%` | `0.944744909` |
| C3D8RT B5.54 | `22.72 s` | `29.508592 s` | `23.0055%` | `0.769945242` |

The corresponding Abaqus analysis wall times are `31 s` and `26 s`. Fuelsim
internal totals are `27.340612002 s` and `18.617592879 s`. C3D8T uses 93
nonlinear iterations, 113 residual evaluations, and 26 Jacobian evaluations;
C3D8RT uses 84, 104, and 35. Both complete all twenty increments without a
rejected step.

The optimization keeps exact current-face projection as the ownership test but
uses the contact search tree to reject faces whose current bounding boxes cannot
be closer. It reuses the prior exact owner only as an initial search bound,
precomputes shared primary edges, evaluates residual contact projection in
ordinary double precision, compresses each active averaged constraint to its
current node support, and permits state-dependent matrix nonzeros only for this
finite-sliding averaged-contact path. Fixed-sparsity problems retain strict
PETSc preallocation and new-nonzero errors.

The complete field comparisons remain separate from timing. For B5.53, contact-
pressure relative L2, relative absolute-peak, and maximum pointwise-relative
errors are `0.0594374%`, `0.0436218%`, and `0.137470%`. For B5.54 they are
`0.0496178%`, `0.0520941%`, and `0.0686526%`. All aggregate field metrics pass
`0.5%`; radial-displacement pointwise errors are `2.90435%` and `2.93536%`
under the existing fixed-path `4%` small-reference qualification. These are
single end-to-end measurements across native Windows and WSL, not medians or
pure kernel ratios.

The tracked Fuelsim inputs are
`verification/fuelsim/transient_integrated_c3d8t_sts.fsi` and
`verification/fuelsim/transient_integrated_c3d8rt_sts.fsi`. The Abaqus decks,
raw reference fields, comparison summaries, and external timing records use
the `b553_m58_c3d8t_sts_*` and `b554_m58_c3d8rt_sts_*` prefixes under
`verification/abaqus`.

## M5.8 C3D20T Abaqus timing and field comparison

This B5.48 node-to-surface result is historical evidence. The tracked input and
timing records are retained, but current Fuelsim builds reject C3D20T
`node_to_surface` mechanical contact during problem construction; current
C3D20T qualification and timing work uses the surface-to-surface B5.51 and
B5.56 paths below.

B5.48 upgrades every element of the tracked M5.8 mesh from HEX8 to HEX20
without changing the 1,152-element partition, cylindrical geometry, two
material regions, contact surfaces, material coefficients, loads, or twenty
fixed `0.05 s` increments. Unique edge midpoints increase the mesh from 1,617
corner nodes to 5,969 displacement nodes. Both Fuelsim and Abaqus activate
temperature only at the 1,617 corners, so each model has 19,524 coupled degrees
of freedom. The conversion keeps Fuelsim's local HEX20 ordering in Exodus and
explicitly maps its top-edge and vertical-edge midpoint groups to the Abaqus
C3D20T ordering in the generated include.

Both solvers use one process and a direct solver. Fuelsim is pinned to CPU 0,
and OpenMP, OpenBLAS, MKL, and NumExpr are each limited to one thread. Abaqus
R2018x uses `cpus=1`. Both complete all twenty increments without a cutback.
Fuelsim takes 86 nonlinear iterations and Abaqus takes 71 total iterations.
The controlled Fuelsim external-wall trials are `295.89 s`, `293.66 s`, and
`291.32 s`; their median is `293.66 s`. The Abaqus external wall time is
`193.928138 s`, with `181.50 s` total CPU time and `190 s` reported wall time
in its job summary. The observed Fuelsim-to-Abaqus external-wall ratio is
`1.514272261`: Fuelsim is `51.427%` slower relative to Abaqus, while Abaqus's
wall time is `33.962%` lower relative to Fuelsim.

This is a same-mesh, same-input-physics end-to-end timing comparison across
Windows and WSL, not a pure element-kernel ratio. The final fields have now
been extracted and compared, and B5.48 does not pass full-field qualification.
Corner-temperature errors remain below `0.025%`, but radial-displacement,
contact-pressure, equivalent-stress, equivalent-plastic-strain, and
equivalent-creep-strain relative L2 errors are `21.9330%`, `29.4997%`,
`37.1520%`, `37.9305%`, and `82.3086%`. The nominal-input runtime ratio must
therefore not be presented as a qualified same-result performance comparison.
The tracked `fuelsim_m58_integrated_hex20_results` target writes one final
Exodus result without changing the timed production input, and
`compare_b548.py` reproduces all three error metrics without a denominator
floor. No production formula or coefficient is changed from this failure.

## B5.51 C3D20T surface-to-surface timing after heat and pressure correction

B5.51 replaces the node-to-surface mechanical interface with frictionless
finite-sliding surface-to-surface contact and constrains all quadratic top and
bottom face displacement nodes in both solvers. It also raises the fuel heat
source to `2e8 W/m^3`. Correct bilinear interpolation of the four C3D20T thermal
`HFL` values on the quadratic current surface gives final contact heat rates of
`14.308063 W` in Fuelsim and `14.342998 W` in Abaqus. This is a
19,524-degree-of-freedom, twenty-step manual benchmark.

On CPU 0 with every listed numerical library restricted to one thread,
Fuelsim takes `317.07 s` externally and `300.873 s` internally. Abaqus R2018x
with `cpus=1` takes `608.050095 s` externally and reports `604 s` analysis wall
time. Fuelsim is `1.918` times as fast by external wall time and `2.007` times as
fast when its internal total is compared with Abaqus analysis time in this
single observation. Fuelsim uses 86 nonlinear iterations, 106 residual
evaluations, 38 Jacobian evaluations, and one PETSc workspace setup; its peak
resident memory is `1,363,652 kB`.

The corrected source-time treatment gives temperature, radial-displacement,
equivalent-stress, equivalent-plastic-strain, and equivalent-creep-strain
relative L2 errors of `0.00848164%`, `0.00336854%`, `0.00683641%`,
`0.00696828%`, and `0.0108461%`. Resolving each curved secondary face's normal
gap before shared-node averaging reduces recovered contact-pressure L2, peak,
and pointwise errors to `0.00209098%`, `0.00298428%`, and `0.00411733%`.
The recovered-pressure surface-integral error is `0.000136010%`; direct radial
projection of native nodal normal-force vectors differs by `0.0000557156%`.
The contact heat-rate difference is `0.243566%`. All matched nonzero physical
fields and contact integrals are below `0.5%`, so B5.51 is qualified. Near-zero
tangential displacement percentages remain diagnostic and use no denominator
floor.

## B5.56 C3D20T surface-to-surface friction timing

B5.56 keeps the complete B5.51 workload and adds Coulomb friction with
`mu=0.002`. Both solvers use finite-sliding surface-to-surface contact and zero
frictional heat conversion. Fuelsim uses one current-geometry constraint and
friction history per unique secondary quadratic node; curved-face samples form
the local normal, tangent, scalar gap, and displacement gradient before their
node-centered average. This reduces the large case from 141,568 contact
contributions on the old nine-points-per-face path to 2,592 while retaining the
full finite-sliding candidate sparsity.

The MUMPS/PORD direct-solver path now keeps that complete assembly matrix and
builds a separate, numerically equal factor matrix. Only exact off-diagonal
zeros are removed, and each Jacobian update rebuilds symbolic factorization.
PORD is the automatic MUMPS ordering because it preserves the tested bitwise
checkpoint replay when repeated symbolic analyses are necessary. Explicit
SCOTCH or other PETSc orderings retain the full-matrix path. No contact law,
material parameter, time step, nonlinear tolerance, or accuracy gate changes.
The diagnostic profile records a reduction in factor setup from `163.54 s` to
`61.457 s`; its output-enabled run is separate from the formal timing samples.

On CPU 0 with one process and every numerical library restricted to one thread,
three production-entry Fuelsim runs take `191.03 s`, `193.43 s`, and `194.99 s`
externally; the median is `193.43 s`. The alternating CPU-0-affinity Abaqus
R2018x `cpus=1` runs take `232.874728 s`, `240.817169 s`, and `239.089524 s`;
their median is `239.089524 s`. Field, history, and restart output are disabled
for formal timing. Fuelsim uses `19.0973%` less external wall time, a `1.23605x`
speedup on this Intel Core i9-13980HX host across Windows and WSL2. Abaqus
analysis wall times are `227`, `236`, and `234 s`; its analysis CPU times are
`218.30`, `223.60`, and `224.00 s`. These remain distinct from external wall
time. The unchanged baseline executable was independently checked at `277.78 s`
with profiling enabled. The formal samples and separate diagnostic profile are
stored alongside the B5.56 reference results.

All 416 secondary contact nodes are projected and sliding. Contact pressure,
gap, normal-force magnitude, tangential-force magnitude, dominant axial shear,
tangential-slip magnitude, temperature, material histories, recovered pressure
and heat integrals pass the tracked `0.5%` gates. The complete tangential
resultant differs by `0.00225030%` and `5.02922e-7 N`, below its explicit
`1e-6 N` absolute gate. The raw Cartesian tangential nodal-force
vector and Abaqus smoothed `CSHEAR` integral remain explicitly reported
diagnostics because they compare unstable near-zero curved-basis components or
different recovered quantities. The full commands and all error rows are in
`verification/abaqus/README.md` and the B5.56 comparison artifact.

## B6.0 C3D8RT fuel-plate thermo-mechanical bending

B6.0 uses a 100 mm × 6 mm × 1 mm plate: 80 central `meat` elements are
continuously enclosed by 1,120 `clad` elements. The structured mesh has 50, 6,
and 4 elements along x, y, and z; the fuel core occupies x=`0.04--0.06 m`,
y=`0.001--0.005 m`, z=`0.00025--0.00075 m`. Both blocks use C3D8RT with the
same four-field degrees of freedom at shared interface nodes. The fuel block
generates `2e8 W/m^3`; the front and back faces are prescribed at `700` and
`600 K`, producing a strong through-thickness temperature gradient and a
measurable free-edge bending displacement. The left face is clamped in all three
displacement components. The path has 1,785 nodes, 1,200 elements, 7,140 degrees
of freedom, and ten fixed one-second increments. The fuel isotropic
thermal-expansion coefficient is `1e-4 K^-1`; the cladding coefficient remains
`5e-6 K^-1`.

Run the Fuelsim side with one process and one thread per numerical library:

```text
./build/fuelsim_b60_fuel_plate_c3d8rt_benchmark \
  verification/fuelsim/transient_b60_fuel_plate_c3d8rt_bending.fsi \
  /tmp/b60_fuelsim_nodal.csv /tmp/b60_fuelsim_timing.tsv
```

On the current machine this completed all ten increments with one PETSc workspace
setup and an internal solver time of `4.36 s` in the representative run. The
benchmark passes the input Jacobian-lag setting through to the solver; for this
linear, constant-coefficient case, `jacobian_lag = 10` reduces Jacobian
evaluations from 125 to 17 without changing the converged field. With the
process fixed to one CPU core and all numerical libraries restricted to one
thread, the external time was `4.45 s` (median of three sequential runs). The final
temperature range is `600--700 K`; the free-right-edge displacement-z range is
`2.4244148e-6 m`, while the largest absolute right-edge z displacement is
`4.474 mm`.

The Abaqus deck uses the same converted mesh include and C3D8RT element,
coupled temperature-displacement transient step, material data, heat source,
clamp, and prescribed front/back temperatures. On a machine with Abaqus R2018x,
run
`verification/abaqus/run_b60.ps1`, then extract and compare:

```text
python3 verification/abaqus/compare_b60.py \
  /tmp/b60_fuelsim_nodal.csv \
  verification/abaqus/b60_fuel_plate_c3d8rt_bending_nodal.csv \
  --fuelsim-timing /tmp/b60_fuelsim_timing.tsv \
  --abaqus-timing verification/abaqus/b60_fuel_plate_c3d8rt_bending_timing.txt \
  --fuelsim-external-seconds 4.45
```

The recorded Abaqus runs used one CPU and full output precision; their external
wall times were `5.394995 s`, `3.387755 s`, and `5.359170 s` (median
`5.359170 s`). Relative L2 errors are `5.86e-10%` for temperature,
`1.85e-6%` for displacement-x, `2.72e-5%` for displacement-y, and
`9.06e-7%` for displacement-z. Relative absolute-peak errors are respectively
`5.65e-9%`, `3.14e-6%`, `3.93e-5%`, and `6.18e-7%`; maximum absolute differences
are `3.95e-8 K`, `3.69e-12 m`, `5.25e-12 m`, and `2.76e-11 m`. The free-edge
bending ranges differ by `5.12e-13 m`. The large pointwise percentages reported
for near-zero displacement references are denominator amplification, not a large
absolute field discrepancy. Using the median external times, Fuelsim/Abaqus is
`0.830`, so Fuelsim is about `1.20x` faster in this cross Windows-and-WSL
measurement.
The reproducible metric summary is stored in
`verification/abaqus/b60_fuel_plate_c3d8rt_bending_comparison.tsv`.

### B6.0 finite-strain diagnostic extension

The same 100 mm mesh and ten one-second increments are also run with
`strain = finite` in both Fuelsim regions and `NLGEOM=YES` in Abaqus. This is a
manual diagnostic comparison, not a qualified production path or a CTest: it
was added to expose the finite-strain C3D8RT mechanics before any accuracy gate
is claimed.

Run Fuelsim with the finite input and run Abaqus with the matching PowerShell
wrapper:

```text
./build/fuelsim_b60_fuel_plate_c3d8rt_benchmark \
  verification/fuelsim/transient_b60_fuel_plate_c3d8rt_finite_bending.fsi \
  /tmp/b60_finite_nodal.csv /tmp/b60_finite_timing.tsv
powershell -ExecutionPolicy Bypass -File verification/abaqus/run_b60_finite.ps1 \
  -SourceDirectory verification/abaqus
python3 verification/abaqus/compare_b60.py \
  /tmp/b60_finite_nodal.csv \
  verification/abaqus/b60_fuel_plate_c3d8rt_finite_bending_nodal.csv \
  --fuelsim-timing /tmp/b60_finite_timing.tsv \
  --abaqus-timing verification/abaqus/b60_fuel_plate_c3d8rt_finite_bending_timing.txt \
  --fuelsim-external-seconds 30.79
```

The controlled three-run medians are `30.79 s` for Fuelsim and `9.425431 s`
for Abaqus, giving a Fuelsim/Abaqus external ratio of `3.26669`; Abaqus is
about `3.27x` faster for this finite-strain path. Temperature relative L2,
relative absolute-peak, and maximum pointwise errors are `0.000517%`,
`0.001409%`, and `0.001461%`, respectively. The corresponding displacement-x,
-y, and -z relative L2 errors are `81.9795%`, `8.64331%`, and `34.6238%`.
The free-right-edge bending ranges are `2.46815e-6 m` (Fuelsim) and
`6.13938e-6 m` (Abaqus), an absolute difference of `3.67123e-6 m`.

These mechanics errors fail the finite-strain comparison gate; the result is
stored as a diagnostic baseline for the next C3D8RT finite-strain correction.
The very large pointwise percentages also include amplification at small but
nonzero displacement references, while the reported absolute differences give
their physical scale. The complete rows are retained in
`verification/abaqus/b60_fuel_plate_c3d8rt_finite_bending_comparison.tsv`.

### B6.0 qualified finite-strain ramped path

The qualified finite-strain comparison retains the same mesh, materials, heat
source, ten one-second increments, C3D8RT elements, and nonlinear-geometry
settings. The only physical-path correction is to prescribe the front-face
temperature as the same linear `600--700 K` history in both solvers. The older
instant-temperature path remains above as a diagnostic because it applies the
full temperature change in the first large-deformation increment and produces
a path-dependent mismatch between the two incremental finite-strain updates.

Run the matched ramped inputs and compare their final nodal fields:

```text
./build/fuelsim_b60_fuel_plate_c3d8rt_benchmark \
  verification/fuelsim/transient_b60_fuel_plate_c3d8rt_finite_ramped_bending.fsi \
  /tmp/b60_finite_ramped_nodal.csv /tmp/b60_finite_ramped_timing.tsv
powershell -ExecutionPolicy Bypass -File verification/abaqus/run_b60_finite_ramped.ps1 \
  -SourceDirectory verification/abaqus
python3 verification/abaqus/compare_b60.py \
  /tmp/b60_finite_ramped_nodal.csv \
  verification/abaqus/b60_fuel_plate_c3d8rt_finite_ramped_bending_nodal.csv \
  --fuelsim-timing /tmp/b60_finite_ramped_timing.tsv \
  --abaqus-timing verification/abaqus/b60_fuel_plate_c3d8rt_finite_ramped_bending_timing.txt \
  --fuelsim-external-seconds 5.51
```

Temperature relative L2, relative absolute-peak, and maximum pointwise errors
are `2.32913e-5%`, `0.000116121%`, and `0.000124927%`. On the free nodes, the
complete displacement-vector errors are `0.00996672%`, `0.0101872%`, and
`0.0444445%`, so all three acceptance metrics are below `0.5%` without a
denominator floor. The X-, Y-, and Z-displacement aggregate errors are also at
most `0.104851%`. Their raw componentwise pointwise percentages remain
diagnostics: the extrema use Abaqus reference components between `7.09e-41 m`
and `5.16e-10 m`. The maximum vector difference is `4.64277e-7 m`, and the
free-right-edge bending ranges differ by `6.19017e-10 m`.

The Fuelsim input uses a linear time predictor, rebuilds the Jacobian at every
Newton iteration, and enables convergence on the same per-field physical
residual tolerances used by the final residual audit. All ten steps converge
without rejection in 21 nonlinear iterations and 21 Jacobian evaluations,
using one PETSc workspace. Abaqus uses 19 nonlinear iterations, 19 equation-
solver passes, and 19 matrix decompositions, so Fuelsim is two iterations, or
`10.5263%`, above Abaqus while remaining at a similar level.

With CPU 0 pinned, MUMPS selected, and all numerical libraries restricted to
one thread, Fuelsim external times were `5.49`, `5.51`, and `5.69 s`, with a
`5.51 s` median. The matching Abaqus `cpus=1` times were `12.649117`,
`7.3038743`, and `7.4664226 s`, with a `7.4664226 s` median; its representative
job summary reports `3.5 s` total CPU time and `4 s` analysis wall time. The
external Fuelsim/Abaqus ratio is `0.737971`, so Fuelsim uses `26.2029%` less
external wall time in this cross-Windows-and-WSL comparison. Complete evidence
is stored in the ramped comparison, timing, and
`b60_fuel_plate_c3d8rt_nonlinear_iteration_diagnosis.tsv` artifacts.

### B6.0 steady-state finite-strain thermoelastic path

The B6.0 benchmark executable also accepts
`verification/fuelsim/steady_b60_fuel_plate_c3d8rt_finite_bending.fsi`. This
input retains the same mesh, C3D8RT elements, fuel heat source, final face
temperatures, material properties, and clamp, but solves one steady load step
without a heat-capacity term. The matching Abaqus input uses one steady-state
coupled temperature-displacement increment with nonlinear geometry enabled.

Centered directional-difference checks at all 11 states visited by the
original zero-displacement Newton path found a worst field-block Jacobian error
of `8.63e-7`. The excess iterations therefore came from applying the full
thermal bending load at a poor initial displacement, not from a missing tangent
term. The optional small-strain solve supplies a close full-load displacement
shape before the unchanged finite-strain equations are solved.

```text
./build/fuelsim_b60_fuel_plate_c3d8rt_benchmark \
  verification/fuelsim/steady_b60_fuel_plate_c3d8rt_finite_bending.fsi \
  /tmp/b60_finite_steady_nodal.csv /tmp/b60_finite_steady_timing.tsv
powershell -ExecutionPolicy Bypass -File verification/abaqus/run_b60_finite_steady.ps1 \
  -SourceDirectory verification/abaqus
python3 verification/abaqus/compare_b60.py \
  /tmp/b60_finite_steady_nodal.csv \
  verification/abaqus/b60_fuel_plate_c3d8rt_finite_steady_bending_nodal.csv
```

Temperature passes relative L2, relative absolute-peak, and maximum pointwise
errors at `2.38862e-5%`, `0.000108914%`, and `0.000111674%`. The corresponding
free-node complete displacement-vector errors are `0.0126146%`, `0.0126317%`,
and `0.0485431%`, with a `5.75351e-7 m` maximum absolute difference. Component
pointwise percentages remain near-zero-reference diagnostics without a
denominator floor. A full-load small-strain predictor takes one iteration, and
the finite-strain corrector takes four more; both reuse one PETSc workspace.
The total is 5 nonlinear iterations and 5 Jacobian evaluations, exactly the
same as Abaqus. The three fixed-CPU external times are `1.27`, `1.30`, and
`1.28 s`; these startup-dominated cross-operating-system timings are diagnostic
only. This comparison qualifies the one-step steady final equilibrium and does
not replace the transient B6.0 path. The baseline, corrected, and Abaqus solver
counts are retained together in
`verification/abaqus/b60_fuel_plate_c3d8rt_nonlinear_iteration_diagnosis.tsv`.

## B6.1 finite-strain inelastic fuel-plate bending

B6.1 retains the 1,785-node, 1,200-element and 7,140-degree-of-freedom B6.0
mesh, heat source, back-face temperature, and clamp. The matched five-step path
ramps the front face from 600 to 800 K and extends the 100 mm plate axially by
0.8 mm. Both finite-strain C3D8RT regions use fully coupled J2 plasticity and
Norton creep. The common verification constants are a 1 MPa initial yield
stress, a 20 GPa hardening modulus, and a Norton rate of `3.5e-4 s^-1` at
100 MPa with exponent 3. Inelastic dissipation is excluded from the heat
equation in both programs. These constants exercise the numerical branches and
are not empirical fuel or cladding models; they remain unchanged in every
performance run.

The registered comparison runs all five increments and checks 8,925 nodal
records, 6,000 material-point records, and five energy records:

```text
ctest --test-dir build \
  -R '^fuelsim_b61_fuel_plate_c3d8rt_finite_inelastic_abaqus_tests$' \
  -j4 --output-on-failure
```

The relative L2, relative absolute-peak, and maximum pointwise-relative errors
are `0.00382850%`, `0.00593139%`, and `0.0678457%` for displacement; `0.0772653%`,
`0.205813%`, and `0.447050%` for reaction force; `0.0194975%`, `0.0690957%`, and
`0.103000%` for stress; `0.00857108%`, `0.000734885%`, and `0.0807925%` for
equivalent plastic strain; and `0.0207731%`, `0.00368103%`, and `0.210174%` for
equivalent creep strain. All use the common 0.5 percent gate without a
denominator floor. All 6,000 material-point rows have nonzero plastic and creep
references. Maximum equivalent plastic and creep strains are `0.00802400` and
`0.00823199`. Abaqus artificial strain energy reaches `0.788155%` of internal
energy.

For the controlled timing runs, CPU 0 was fixed, MUMPS was selected, and every
numerical library was limited to one thread. Fuelsim production-entry external
times were `5.30`, `5.34`, and `5.34 s`, with a `5.34 s` median. Abaqus
`cpus=1` external times were `5.389362`, `5.4808106`, and `5.3975893 s`, with a
`5.3975893 s` median. The Fuelsim-to-Abaqus ratio is `0.989331`, so Fuelsim uses
`1.06694%` less external wall time in this cross-Windows-and-WSL comparison.
The same physical workload through the production entry at pre-change commit
`26dde87`, using its default SCOTCH ordering, has a `5.88 s` median; the current
PORD-selected path is `9.18367%` lower without a material change.
For a stricter code-path isolation, forcing PORD on both commits gives `5.35`
and `5.34 s` medians, a `0.186916%` reduction from the C++ assembly and data
changes alone. Abaqus
uses 14 nonlinear iterations and matrix decompositions; Fuelsim uses 26
nonlinear iterations and Jacobians with one PETSc workspace. External wall
time, not either program's internal timer, is the speed comparison metric.
Complete field, energy, solver, and timing data are retained in
`verification/abaqus/b61_fuel_plate_c3d8rt_finite_inelastic_bending_comparison.tsv`.

The required unrelated solver regression used the same CPU and thread limits.
For the 1,584-degree-of-freedom case, pre-change commit `26dde87` formal samples
were `1.42`, `1.41`, and `1.39 s`, while the current samples were `1.40`, `1.40`,
and `1.38 s`; the medians are `1.41` and `1.40 s`, respectively, with 64
nonlinear iterations on both sides. The current 23,010-degree-of-freedom,
20-step direct MUMPS case completed all 20 steps in `27.81 s` externally and
`27.3085 s` internally, with 62 nonlinear iterations, 82 residual evaluations,
62 Jacobian evaluations, one PETSc workspace, and a final residual norm of
`3.24913e-9`. These runs check for a general solver regression and do not replace
the matched B6.1 Abaqus comparison.

## Manual C3D20T finite-strain plate bending

Three C3D20T cases use the 1,200-element quadratic-displacement and
linear-temperature plate mesh: one steady thermoelastic equilibrium, ten fixed
one-second thermoelastic increments, and five fixed two-second increments with
fully coupled J2 plasticity and Norton creep. They are intentionally manual and
are not registered in CTest.

The finite-strain thermal volume operators match independently identified
Abaqus behavior: 27-point conduction with midpoint gradients and current-volume
measure, 27-point consistent current-volume heat capacity, and current-volume
body source based on the eight temperature corners. The geometric Jacobian uses
closed double-precision derivatives and reuses the mechanical kinematics; no
material coefficient or time step was changed for performance.

For fixed-core, one-thread external wall time, the Fuelsim and Abaqus medians
are `14.67` and `19.555856 s` for steady bending, `57.40` and `75.822728 s` for
the ten-step path, and `52.75` and `55.964890 s` for the plasticity-creep path.
Fuelsim is respectively `24.9841%`, `24.2971%`, and `5.74448%` faster.

The plasticity-creep case passes all three `0.5%` field metrics, including all
32,400 integration-point stress, equivalent plastic-strain, and equivalent
creep-strain values. The thermoelastic cases pass temperature and complete
displacement-vector metrics, while maximum pointwise stress errors remain
`2.47082%` and `15.3961%` only at kilopascal-scale reference stresses; their
relative L2 errors are `3.17327e-5%` and `0.000264014%`. No denominator floor or
extra accuracy treatment is used. These two pointwise diagnostics are qualified
only for the recorded cases with respective `2.5%` and `15.5%` boundaries; all
other field metrics retain the `0.5%` boundary. Complete commands and evidence are in
`verification/abaqus/README.md` and the corresponding comparison and timing
artifacts.
