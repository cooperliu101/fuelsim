[Case]
  version = 3
  problem = steady
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../abaqus/b41_hex8_finite_sliding_mesh.e
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
  [primary]
    block = primary
    material = elastic
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary]
    block = secondary
    material = elastic
    strain = finite
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
      sliding = finite
      penalty = 1e5
      mu = 0.5
      elastic_slip = 1e-4
    []
  []
[]

[BoundaryConditions]
  [primary_temperature]
    type = dirichlet
    boundary = primary_all
    field = temperature
    value = 300
  []
  [secondary_temperature]
    type = dirichlet
    boundary = secondary_all
    field = temperature
    value = 300
  []
  [primary_x]
    type = dirichlet
    boundary = primary_all
    field = displacement_x
    value = 0
  []
  [primary_y]
    type = dirichlet
    boundary = primary_all
    field = displacement_y
    value = 0
  []
  [primary_z]
    type = dirichlet
    boundary = primary_all
    field = displacement_z
    value = 0
  []
  [secondary_x]
    type = dirichlet
    boundary = secondary_all
    field = displacement_x
    value = -0.01
    scale_with_load = true
  []
  [secondary_y]
    type = dirichlet
    boundary = secondary_all
    field = displacement_y
    value = 1
    scale_with_load = true
  []
  [secondary_z]
    type = dirichlet
    boundary = secondary_all
    field = displacement_z
    value = 0.4
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
[]
