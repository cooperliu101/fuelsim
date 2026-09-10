[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../abaqus/h20_24_hex20_nonmatching_contact_mesh.e
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
    element = c3d20t
    block = primary
    material = primary
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary]
    element = c3d20t
    block = secondary
    material = secondary
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]

[TimeFunctions]
  [normal_path]
    type = piecewise_linear
    times = 0 1 2 3 4 5 6 7
    values = 0 -1e-5 -1e-5 -1e-5 -1e-5 -1e-5 -1e-5 -1e-5
  []
  [tangential_path]
    type = piecewise_linear
    times = 0 1 2 3 4 5 6 7
    values = 0 5e-6 1e-5 3e-5 6e-5 1e-5 -5e-5 -4.5e-5
  []
[]

[Contact]
  [interface]
    primary = primary_contact
    secondary = secondary_contact
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
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
    value = 300
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
  end_time = 7
  initial_time_step = 0.1
  minimum_time_step = 0.1
  maximum_time_step = 0.1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 3
  load_ramp_time = 0
[]

[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-11
  step_tolerance = 1e-12
  maximum_iterations = 40
  linear_solver = direct
  direct_factorization = mumps
  field_residual_scaling = true
  linear_relative_tolerance = 1e-11
  maximum_linear_iterations = 400
[]

[Outputs]
  console = true
  csv = transient_hex20_nonmatching_friction_abaqus_summary.csv
  exodus = transient_hex20_nonmatching_friction_abaqus_results.e
  exodus_interval = 10
[]
