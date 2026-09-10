[Case]
  version = 3
  problem = steady
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

[Contact]
  [interface]
    primary = primary_contact
    secondary = secondary_contact
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      penalty = 1e11
      mu = 0
    []
  []
[]

[BoundaryConditions]
  [primary_temperature]
    type = dirichlet
    boundary = primary_z0
    field = temperature
    value = 300
  []
  [secondary_temperature]
    type = dirichlet
    boundary = secondary_z0
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
    value = -1e-5
    scale_with_load = true
  []
  [secondary_y]
    type = dirichlet
    boundary = secondary_y0
    field = displacement_y
    value = 0
  []
  [secondary_z]
    type = dirichlet
    boundary = secondary_z0
    field = displacement_z
    value = 0
  []
[]

[Executioner]
  type = steady
  load_steps = 1
[]

[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-11
  step_tolerance = 1e-12
  maximum_iterations = 30
  linear_solver = direct
  direct_factorization = mumps
  field_residual_scaling = true
  linear_relative_tolerance = 1e-11
  maximum_linear_iterations = 400
[]

[Outputs]
  console = true
  csv = steady_hex20_nonmatching_contact_abaqus_summary.csv
  exodus = steady_hex20_nonmatching_contact_abaqus_results.e
[]
