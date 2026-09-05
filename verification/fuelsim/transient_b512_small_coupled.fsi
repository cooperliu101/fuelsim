[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = ../meshes/b510_hex8_unit_cube.e
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 1
      density = 1
      specific_heat = 1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.3
    []
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 2e8
      hardening_modulus = 2e9
    []
    [creep]
      function = norton
      coefficient = 1e-4
      reference_stress = 1e8
      stress_exponent = 3
    []
  []
[]
[Regions]
  [solid]
    block = solid
    material = solid
    strain = small
    element = c3d8t
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[TimeFunctions]
  [right_displacement]
    type = piecewise_linear
    times = 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
    values = 0 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005 0.0055 0.006
  []
[]
[BoundaryConditions]
  [temperature]
    type = dirichlet
    boundary = left
    field = temperature
    value = 600
  []
  [fix_x]
    type = dirichlet
    boundary = left
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = y0
    field = displacement_y
    value = 0
  []
  [fix_z]
    type = dirichlet
    boundary = z0
    field = displacement_z
    value = 0
  []
  [pull_x]
    type = dirichlet
    boundary = right
    field = displacement_x
    value = 1
    function = right_displacement
  []
[]
[Executioner]
  type = transient
  end_time = 1
  initial_time_step = 0.1
  minimum_time_step = 0.1
  maximum_time_step = 0.1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-12
  step_tolerance = 1e-12
  maximum_iterations = 40
  linear_solver = direct
  preconditioner = lu
  line_search = basic
  backtracking_fallback = true
  field_residual_scaling = true
  residual_reduction_tolerance = 1e-11
  temperature_residual_absolute_tolerance = 1e-10
  mechanical_residual_absolute_tolerance = 1e-2
[]
[Outputs]
  console = true
  csv = transient_b512_small_coupled_summary.csv
  exodus = transient_b512_small_coupled_results.e
[]
