[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/b114_cax8t_recovery.e
[]
[TimeFunctions]
  [r0]
    type = piecewise_linear
    times = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
    values = 0 0.000002 0.000006 0.000002 0.000002 0.000002 0.000002 0.000003 0 0 0 0 0 0.000002 0.000002
  []
  [z0]
    type = piecewise_linear
    times = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
    values = 0 0 0 0 0 0 0 0.000001 0.000001 0.000001 0.000001 0.000001 0 1e-8 0.000002
  []
  [r1]
    type = piecewise_linear
    times = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
    values = 0 0.000002 0.000002 0.000006 0.000002 0.000002 0.000002 0 0.000003 0 0 0 0 0.0000024999999999999998 0.0000024999999999999998
  []
  [z1]
    type = piecewise_linear
    times = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
    values = 0 0 0 0 0 0 0 0.000001 0.000001 0.000001 0.000001 0.000001 0 -2e-8 0.000002
  []
  [r2]
    type = piecewise_linear
    times = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
    values = 0 0.000002 0.000002 0.000002 0.000006 0.000002 0.000002 0 0 0.000003 0 0 0 0.000003 0.000003
  []
  [z2]
    type = piecewise_linear
    times = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
    values = 0 0 0 0 0 0 0 0.000001 0.000001 0.000001 0.000001 0.000001 0 3.0000000000000004e-8 0.000002
  []
  [r3]
    type = piecewise_linear
    times = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
    values = 0 0.000002 0.000002 0.000002 0.000002 0.000006 0.000002 0 0 0 0.000003 0 0 0.0000035 0.0000035
  []
  [z3]
    type = piecewise_linear
    times = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
    values = 0 0 0 0 0 0 0 0.000001 0.000001 0.000001 0.000001 0.000001 0 -4e-8 0.000002
  []
  [r4]
    type = piecewise_linear
    times = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
    values = 0 0.000002 0.000002 0.000002 0.000002 0.000002 0.000006 0 0 0 0 0.000003 0 0.000004 0.000004
  []
  [z4]
    type = piecewise_linear
    times = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
    values = 0 0 0 0 0 0 0 0.000001 0.000001 0.000001 0.000001 0.000001 0 5e-8 0.000002
  []
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.3
    []
  []
[]
[Regions]
  [inner]
    block = inner
    element = cax8rt
    material = solid
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
  [outer]
    block = outer
    element = cax8rt
    material = solid
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[Contact]
  [interface]
    primary = outer_left
    secondary = inner_right
    [mechanical]
      formulation = penalty
      discretization = node_to_surface
      sliding = finite
      penalty = 1e13
      mu = 0.2
      slip_tolerance = 0.001
    []
  []
[]
[BoundaryConditions]
  [temperature]
    type = dirichlet
    boundary = all
    field = temperature
    value = 600
  []
  [fixed_r]
    type = dirichlet
    boundary = fixed
    field = radial_displacement
    value = 0
  []
  [fixed_z]
    type = dirichlet
    boundary = fixed
    field = axial_displacement
    value = 0
  []
  [r0]
    type = dirichlet
    boundary = c0
    field = radial_displacement
    value = 1
    function = r0
  []
  [z0]
    type = dirichlet
    boundary = c0
    field = axial_displacement
    value = 1
    function = z0
  []
  [r1]
    type = dirichlet
    boundary = c1
    field = radial_displacement
    value = 1
    function = r1
  []
  [z1]
    type = dirichlet
    boundary = c1
    field = axial_displacement
    value = 1
    function = z1
  []
  [r2]
    type = dirichlet
    boundary = c2
    field = radial_displacement
    value = 1
    function = r2
  []
  [z2]
    type = dirichlet
    boundary = c2
    field = axial_displacement
    value = 1
    function = z2
  []
  [r3]
    type = dirichlet
    boundary = c3
    field = radial_displacement
    value = 1
    function = r3
  []
  [z3]
    type = dirichlet
    boundary = c3
    field = axial_displacement
    value = 1
    function = z3
  []
  [r4]
    type = dirichlet
    boundary = c4
    field = radial_displacement
    value = 1
    function = r4
  []
  [z4]
    type = dirichlet
    boundary = c4
    field = axial_displacement
    value = 1
    function = z4
  []
[]
[Executioner]
  type = transient
  end_time = 1.4
  initial_time_step = 0.1
  minimum_time_step = 0.1
  maximum_time_step = 0.1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-10
  maximum_iterations = 60
[]
[Outputs]
  console = true
  csv = transient_b1214_cax8rt_recovery_summary.csv
  exodus = transient_b1214_cax8rt_recovery_results.e
[]
