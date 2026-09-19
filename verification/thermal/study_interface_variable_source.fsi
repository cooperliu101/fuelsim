[Case]
  version = 3
  physics = thermal
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = study_dc3d8_interface.e
[]
[TimeFunctions]
  [power]
    type = piecewise_linear
    times = 0 1 2
    values = 1 2 0.5
  []
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
    heat_source_function = power
    heat_source_time_evaluation = interval_average
  []
  [cold]
    block = cold
    element = dc3d8
    material = second
    initial_temperature = 300
    volumetric_heat_source = 2000000
    heat_source_function = power
    heat_source_time_evaluation = interval_average
  []
[]
[Contact]
  [gap]
    primary = primary
    secondary = secondary
    [thermal]
      law = affine
      conductance = 1000
      discretization = surface_to_surface
    []
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
  exodus = study_interface_variable_source_results.e
  history = study_interface_variable_source_history.csv
[]
