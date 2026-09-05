[Case]
  version = 3
  problem = steady
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../abaqus/h20_29_hex20_sts_graded_mesh.e
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
      young_modulus = 2000000000
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
      young_modulus = 650000000
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
      penalty = 100000000000
      mu = 0
    []
  []
[]

[BoundaryConditions]
  [primary_temperature]
    type = dirichlet
    boundary = primary_outer
    field = temperature
    value = 300
  []
  [secondary_temperature]
    type = dirichlet
    boundary = secondary_outer
    field = temperature
    value = 300
  []
  [primary_displacement_0]
    type = dirichlet
    boundary = primary_outer
    field = displacement_x
    value = 0
  []
  [primary_displacement_1]
    type = dirichlet
    boundary = primary_outer
    field = displacement_y
    value = 0
  []
  [primary_displacement_2]
    type = dirichlet
    boundary = primary_outer
    field = displacement_z
    value = 0
  []
  [secondary_displacement_0_0]
    type = dirichlet
    boundary = secondary_outer
    field = displacement_x
    value = -1e-05
    scale_with_load = true
  []
  [secondary_displacement_0_1]
    type = dirichlet
    boundary = secondary_outer
    field = displacement_y
    value = -0
    scale_with_load = true
  []
  [secondary_displacement_0_2]
    type = dirichlet
    boundary = secondary_outer
    field = displacement_z
    value = 0
    scale_with_load = true
  []
  [secondary_top_traction]
    type = traction
    boundary = secondary_top
    field = displacement_x
    value = -1000
    configuration = reference
    scale_with_load = true
  []
[]

[Executioner]
  type = steady
  load_steps = 4
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
  csv = steady_hex20_sts_graded_abaqus_summary.csv
  exodus = steady_hex20_sts_graded_abaqus_results.e
[]
