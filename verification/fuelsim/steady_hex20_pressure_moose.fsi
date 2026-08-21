[Case]
  version = 3
  problem = steady
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../moose/h20_03_pressure.e
[]

[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 6000
      specific_heat = 1000
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.25
    []
  []
[]

[Regions]
  [solid]
    block_id = 0
    material = solid
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]

[BoundaryConditions]
  [temperature]
    type = dirichlet
    boundary = left
    field = temperature
    value = 300
  []
  [fix_x]
    type = dirichlet
    boundary = left
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = bottom
    field = displacement_y
    value = 0
  []
  [fix_z]
    type = dirichlet
    boundary = back
    field = displacement_z
    value = 0
  []
  [pressure]
    type = pressure
    boundary = right
    value = 1000000
    configuration = reference
  []
[]

[Executioner]
  type = steady
  load_steps = 1
  cutback_factor = 0.5
  maximum_cutbacks = 3
  minimum_load_increment = 1e-6
[]

[Solver]
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-12
  step_tolerance = 1e-12
  maximum_iterations = 20
[]

[Outputs]
  console = false
[]
