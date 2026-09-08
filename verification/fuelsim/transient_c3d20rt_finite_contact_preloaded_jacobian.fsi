[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../meshes/c3d20rt_finite_contact.e
[]

[Materials]
  [primary]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1
      specific_heat = 1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0
    []
  []
  [secondary]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1
      specific_heat = 1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0
    []
  []
[]

[Regions]
  [primary]
    block = primary
    element = c3d20rt
    material = primary
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary]
    block = secondary
    element = c3d20rt
    material = secondary
    strain = finite
    initial_temperature = 400
    volumetric_heat_source = 0
  []
[]

[TimeFunctions]
  [normal_path]
    type = piecewise_linear
    times = 0 1 2 3 4 5 6 7
    values = 0 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02
  []
  [tangential_path]
    type = piecewise_linear
    times = 0 1 2 3 4 5 6 7
    values = 0 0.005 0.01 0.03 0.06 0.01 -0.05 -0.045
  []
[]

[Contact]
  [interface]
    primary = primary_contact
    secondary = secondary_contact
    [thermal]
      gap_conductivity = 0.001
      minimum_gap = 1e-5
    []
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      sliding = finite
      penalty = 1e11
      mu = 0.3
      slip_tolerance = 1e-5
    []
  []
[]

[BoundaryConditions]
  [primary_temperature]
    type = dirichlet
    boundary = primary_x0
    field = temperature
    value = 300
  []
  [secondary_temperature]
    type = dirichlet
    boundary = secondary_x2
    field = temperature
    value = 400
  []
  [primary_x]
    type = dirichlet
    boundary = primary_x0
    field = displacement_x
    value = 0
  []
  [primary_y]
    type = dirichlet
    boundary = primary_x0
    field = displacement_y
    value = 0
  []
  [primary_z]
    type = dirichlet
    boundary = primary_x0
    field = displacement_z
    value = 0
  []
  [secondary_x]
    type = dirichlet
    boundary = secondary_x2
    field = displacement_x
    value = 1
    function = normal_path
  []
  [secondary_y]
    type = dirichlet
    boundary = secondary_x2
    field = displacement_y
    value = 1
    function = tangential_path
  []
  [secondary_z]
    type = dirichlet
    boundary = secondary_z0
    field = displacement_z
    value = 0
  []
[]

[Executioner]
  type = transient
  end_time = 0.05
  restart = c3d20rt_finite_contact_first_increment.checkpoint
  initial_time_step = 0.025
  minimum_time_step = 0.025
  maximum_time_step = 0.025
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 3
  load_ramp_time = 0
[]

[Solver]
  absolute_tolerance = 1e-12
  relative_tolerance = 1e-13
  step_tolerance = 1e-12
  maximum_iterations = 40
  line_search = backtracking
  linear_solver = direct
  direct_factorization = mumps
  temperature_residual_scale = 1000
  mechanical_residual_scale = 20000000
  linear_relative_tolerance = 1e-11
  maximum_linear_iterations = 400
[]

[Outputs]
  console = true
  csv = transient_c3d20rt_finite_contact_preloaded_jacobian_summary.csv
  exodus = transient_c3d20rt_finite_contact_preloaded_jacobian_results.e
  exodus_interval = 1
[]
