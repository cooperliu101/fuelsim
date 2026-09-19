[Case]
  version = 3
  physics = thermal
  problem = steady
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = coupled.e
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
  [clad]
    [thermal]
      function = constant_thermophysical
      conductivity = 16
      density = 6500
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
    pellet_response = full
  []
  [cladding]
    block = cladding
    element = dc3d8
    material = clad
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[Contact]
  [gap]
    primary = clad_inner
    secondary = pellet_side
    [thermal]
      law = affine
      conductance = 5000
      discretization = surface_to_surface
    []
  []
[]
[BoundaryConditions]
  [coolant]
    type = convection
    boundary = clad_outer
    heat_transfer_coefficient = 10000
    ambient_temperature = 550
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
  exodus = coupled_full_results.e
  csv = coupled_full_history.csv
[]
