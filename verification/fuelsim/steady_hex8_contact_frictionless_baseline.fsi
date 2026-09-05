[Case]
  version = 3
  problem = steady
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../moose/b33_hex8_contact_mesh.e
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
      poisson_ratio = 0.25
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
      poisson_ratio = 0.25
    []
  []
[]

[Regions]
  [primary]
    block = primary
    material = primary
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary]
    block = secondary
    material = secondary
    strain = small
    initial_temperature = 400
    volumetric_heat_source = 0
  []
[]

[Contact]
  [interface]
    primary = primary_right
    secondary = secondary_left
    [thermal]
      gap_conductivity = 0.2
      minimum_gap = 1e-5
    []
    [mechanical]
      formulation = penalty
      penalty = 1e13
      mu = 0
    []
  []
[]

[BoundaryConditions]
  [primary_temperature]
    type = dirichlet
    boundary = primary_left
    field = temperature
    value = 300
  []
  [primary_x]
    type = dirichlet
    boundary = primary_left
    field = displacement_x
    value = 0
  []
  [primary_y]
    type = dirichlet
    boundary = primary_left
    field = displacement_y
    value = 0
  []
  [primary_z]
    type = dirichlet
    boundary = primary_left
    field = displacement_z
    value = 0
  []
  [secondary_temperature]
    type = dirichlet
    boundary = secondary_right
    field = temperature
    value = 400
  []
  [secondary_x]
    type = dirichlet
    boundary = secondary_right
    field = displacement_x
    value = -2e-4
    scale_with_load = true
  []
  [secondary_y]
    type = dirichlet
    boundary = secondary_right
    field = displacement_y
    value = 5e-5
    scale_with_load = true
  []
  [secondary_z]
    type = dirichlet
    boundary = secondary_right
    field = displacement_z
    value = 0
    scale_with_load = true
  []
[]

[Executioner]
  type = steady
  load_steps = 10
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
  csv = steady_hex8_contact_frictionless_baseline_summary.csv
  exodus = steady_hex8_contact_frictionless_baseline_results.e
[]
