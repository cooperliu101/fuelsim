[Case]
  version = 3
  physics = thermal
  problem = steady
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = pellet.e
[]
[Materials]
  [fuel]
    [thermal]
      function = constant_thermophysical
      conductivity = 3
      density = 10000
      specific_heat = 300
    []
  []
[]
[Regions]
  [pellet]
    block = pellet
    element = dc3d8
    material = fuel
    initial_temperature = 600
    volumetric_heat_source = 200000000
    pellet_response = surrogate
    pellet_model = ../models/pellet_mlp.txt
  []
[]
[BoundaryConditions]
  [surface]
    type = dirichlet
    boundary = pellet_boundary
    field = temperature
    value = 600
  []
[]
[Executioner]
  type = steady
  load_steps = 1
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 30
[]
[Outputs]
  console = true
  exodus = pellet_surrogate_results.e
[]
