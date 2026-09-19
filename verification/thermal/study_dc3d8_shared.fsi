[Case]
  version = 3
  physics = thermal
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = study_dc3d8_shared.e
[]
[Materials]
  [first]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
    []
  []
  [second]
    [thermal]
      function = constant_thermophysical
      conductivity = 20
      density = 2000
      specific_heat = 200
    []
  []
[]
[Regions]
  [hot]
    block = hot
    element = dc3d8
    material = first
    initial_temperature = 300
    volumetric_heat_source = 1000000
  []
  [cold]
    block = cold
    element = dc3d8
    material = second
    initial_temperature = 300
    volumetric_heat_source = 2000000
  []
[]
[BoundaryConditions]
  [hot]
    type = dirichlet
    boundary = left
    field = temperature
    value = 400
  []
  [cold]
    type = convection
    boundary = right
    heat_transfer_coefficient = 100
    ambient_temperature = 300
  []
[]
[Executioner]
  type = transient
  end_time = 2
  initial_time_step = 0.1
  minimum_time_step = 0.1
  maximum_time_step = 0.1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 40
[]
[Outputs]
  console = true
  exodus = study_dc3d8_shared_results.e
  history = study_dc3d8_shared_history.csv
[]
