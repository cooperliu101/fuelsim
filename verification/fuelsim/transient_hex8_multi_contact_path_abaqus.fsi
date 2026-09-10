[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../abaqus/b48_hex8_multi_contact_mesh.e
[]

[Materials]
  [elastic]
    [thermal]
      function = constant_thermophysical
      conductivity = 1
      density = 1
      specific_heat = 1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.25
    []
  []
[]

[Regions]
  [primary_a]
    element = c3d8t
    block = primary_a
    material = elastic
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary_a]
    element = c3d8t
    block = secondary_a
    material = elastic
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [primary_b]
    element = c3d8t
    block = primary_b
    material = elastic
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary_b]
    element = c3d8t
    block = secondary_b
    material = elastic
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]

[TimeFunctions]
  [pair_a_x]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 -0.01 -0.01 -0.01 -0.01
  []
  [pair_a_y]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.012 0.36 0.54 0.72
  []
  [pair_a_z]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.014 0.42 0.63 0.84
  []
  [pair_b_x]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 -0.01 -0.01 -0.01 -0.01
  []
  [pair_b_y]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 -0.012 -0.36 -0.54 -0.72
  []
  [pair_b_z]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.014 0.42 0.63 0.84
  []
[]

[Contact]
  [pair_a]
    primary = primary_a_contact
    secondary = secondary_a_contact
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      sliding = finite
      penalty = 1e5
      mu = 0.3
      slip_tolerance = 0.025
    []
  []
  [pair_b]
    primary = primary_b_contact
    secondary = secondary_b_contact
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      sliding = finite
      penalty = 1e5
      mu = 0.5
      slip_tolerance = 0.025
    []
  []
[]

[BoundaryConditions]
  [primary_a_temperature]
    type = dirichlet
    boundary = primary_a_all
    field = temperature
    value = 300
  []
  [secondary_a_temperature]
    type = dirichlet
    boundary = secondary_a_all
    field = temperature
    value = 300
  []
  [primary_b_temperature]
    type = dirichlet
    boundary = primary_b_all
    field = temperature
    value = 300
  []
  [secondary_b_temperature]
    type = dirichlet
    boundary = secondary_b_all
    field = temperature
    value = 300
  []
  [primary_a_x]
    type = dirichlet
    boundary = primary_a_all
    field = displacement_x
    value = 0
  []
  [primary_a_y]
    type = dirichlet
    boundary = primary_a_all
    field = displacement_y
    value = 0
  []
  [primary_a_z]
    type = dirichlet
    boundary = primary_a_all
    field = displacement_z
    value = 0
  []
  [primary_b_x]
    type = dirichlet
    boundary = primary_b_all
    field = displacement_x
    value = 0
  []
  [primary_b_y]
    type = dirichlet
    boundary = primary_b_all
    field = displacement_y
    value = 0
  []
  [primary_b_z]
    type = dirichlet
    boundary = primary_b_all
    field = displacement_z
    value = 0
  []
  [secondary_a_x]
    type = dirichlet
    boundary = secondary_a_all
    field = displacement_x
    value = 1
    function = pair_a_x
  []
  [secondary_a_y]
    type = dirichlet
    boundary = secondary_a_all
    field = displacement_y
    value = 1
    function = pair_a_y
  []
  [secondary_a_z]
    type = dirichlet
    boundary = secondary_a_all
    field = displacement_z
    value = 1
    function = pair_a_z
  []
  [secondary_b_x]
    type = dirichlet
    boundary = secondary_b_all
    field = displacement_x
    value = 1
    function = pair_b_x
  []
  [secondary_b_y]
    type = dirichlet
    boundary = secondary_b_all
    field = displacement_y
    value = 1
    function = pair_b_y
  []
  [secondary_b_z]
    type = dirichlet
    boundary = secondary_b_all
    field = displacement_z
    value = 1
    function = pair_b_z
  []
[]

[Executioner]
  type = transient
  end_time = 4
  initial_time_step = 1
  minimum_time_step = 0.125
  maximum_time_step = 1
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
  csv = transient_hex8_multi_contact_path_abaqus_summary.csv
  exodus = transient_hex8_multi_contact_path_abaqus_results.e
  checkpoint = b48_full.checkpoint
  checkpoint_interval = 1
[]
