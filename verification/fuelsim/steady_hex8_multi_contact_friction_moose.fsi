[Case]
  version = 3
  problem = steady
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../moose/b37_hex8_multi_contact_mesh.e
[]

[Materials]
  [solid]
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
  [primary_a]
    element = c3d8t
    block = primary_a
    material = solid
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary_a]
    element = c3d8t
    block = secondary_a
    material = solid
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [primary_b]
    element = c3d8t
    block = primary_b
    material = solid
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary_b]
    element = c3d8t
    block = secondary_b
    material = solid
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]

[Contact]
  [pair_a]
    primary = primary_a_right
    secondary = secondary_a_left
    [mechanical]
      formulation = penalty
      penalty = 1e12
      mu = 0.001
    []
  []
  [pair_b]
    primary = primary_b_right
    secondary = secondary_b_left
    [mechanical]
      formulation = penalty
      penalty = 1e12
      mu = 0.001
    []
  []
[]

[BoundaryConditions]
  [primary_a_temperature_left]
    type = dirichlet
    boundary = primary_a_left
    field = temperature
    value = 300
  []
  [primary_a_temperature_right]
    type = dirichlet
    boundary = primary_a_right
    field = temperature
    value = 300
  []
  [primary_a_x]
    type = dirichlet
    boundary = primary_a_left
    field = displacement_x
    value = 0
  []
  [primary_a_y]
    type = dirichlet
    boundary = primary_a_left
    field = displacement_y
    value = 0
  []
  [primary_a_z]
    type = dirichlet
    boundary = primary_a_left
    field = displacement_z
    value = 0
  []
  [secondary_a_temperature_left]
    type = dirichlet
    boundary = secondary_a_left
    field = temperature
    value = 300
  []
  [secondary_a_temperature_right]
    type = dirichlet
    boundary = secondary_a_right
    field = temperature
    value = 300
  []
  [secondary_a_x]
    type = dirichlet
    boundary = secondary_a_right
    field = displacement_x
    value = -2e-6
    scale_with_load = true
  []
  [secondary_a_y]
    type = dirichlet
    boundary = secondary_a_right
    field = displacement_y
    value = 2e-6
    scale_with_load = true
  []
  [secondary_a_z]
    type = dirichlet
    boundary = secondary_a_right
    field = displacement_z
    value = 0
    scale_with_load = true
  []
  [primary_b_temperature_left]
    type = dirichlet
    boundary = primary_b_left
    field = temperature
    value = 300
  []
  [primary_b_temperature_right]
    type = dirichlet
    boundary = primary_b_right
    field = temperature
    value = 300
  []
  [primary_b_x]
    type = dirichlet
    boundary = primary_b_left
    field = displacement_x
    value = 0
  []
  [primary_b_y]
    type = dirichlet
    boundary = primary_b_left
    field = displacement_y
    value = 0
  []
  [primary_b_z]
    type = dirichlet
    boundary = primary_b_left
    field = displacement_z
    value = 0
  []
  [secondary_b_temperature_left]
    type = dirichlet
    boundary = secondary_b_left
    field = temperature
    value = 300
  []
  [secondary_b_temperature_right]
    type = dirichlet
    boundary = secondary_b_right
    field = temperature
    value = 300
  []
  [secondary_b_x]
    type = dirichlet
    boundary = secondary_b_right
    field = displacement_x
    value = -2e-6
    scale_with_load = true
  []
  [secondary_b_y]
    type = dirichlet
    boundary = secondary_b_right
    field = displacement_y
    value = 2e-6
    scale_with_load = true
  []
  [secondary_b_z]
    type = dirichlet
    boundary = secondary_b_right
    field = displacement_z
    value = 0
    scale_with_load = true
  []
[]

[Executioner]
  type = steady
  load_steps = 1
[]

[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 50
  temperature_residual_scale = 1e-2
  mechanical_residual_scale = 1e-3
[]

[Outputs]
  console = true
  csv = steady_hex8_multi_contact_friction_moose_summary.csv
  exodus = steady_hex8_multi_contact_friction_moose_results.e
[]
