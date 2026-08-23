[Case]
  version = 3
  problem = steady
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../abaqus/h20_30_hex20_sts_quadratic_cylinder_mesh.e
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
    boundary = primary_inner
    field = temperature
    value = 300
  []
  [secondary_temperature]
    type = dirichlet
    boundary = secondary_outer
    field = temperature
    value = 300
  []
  [primary_fix_x]
    type = dirichlet
    boundary = primary_inner
    field = displacement_x
    value = 0
  []
  [primary_fix_y]
    type = dirichlet
    boundary = primary_inner
    field = displacement_y
    value = 0
  []
  [primary_fix_z]
    type = dirichlet
    boundary = primary_back
    field = displacement_z
    value = 0
  []
  [secondary_fix_z]
    type = dirichlet
    boundary = secondary_back
    field = displacement_z
    value = 0
  []
  [secondary_symmetry_lower]
    type = dirichlet
    boundary = secondary_theta_lower
    field = displacement_y
    value = 0
  []
  [secondary_symmetry_upper]
    type = dirichlet
    boundary = secondary_theta_upper
    field = displacement_x
    value = 0
  []
  [outer_pressure]
    type = pressure
    boundary = secondary_outer
    value = 10000
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
[]
