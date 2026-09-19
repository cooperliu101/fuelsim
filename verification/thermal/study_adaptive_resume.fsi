[Case]
  version = 3
  physics = thermal
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = dc3d20.e
[]
[TimeFunctions]
  [power]
    type = piecewise_linear
    times = 0 1 2
    values = 1 2 0.5
  []
[]
[Materials]
  [solid]
    [thermal]
      function = linear_temperature_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
      reference_temperature = 300
      conductivity_temperature_coefficient = 0
      specific_heat_temperature_coefficient = 1
    []
  []
[]
[Regions]
  [solid]
    block = solid
    element = dc3d20
    material = solid
    initial_temperature = 300
    volumetric_heat_source = 1000000
    heat_source_function = power
  []
[]
[BoundaryConditions]
[]
[Executioner]
  type = transient
  end_time = 2
  initial_time_step = 1
  minimum_time_step = 0.0001
  maximum_time_step = 1
  growth_factor = 2
  cutback_factor = 0.5
  maximum_cutbacks = 20
  load_ramp_time = 0
  adaptive_algorithm = step_doubling
  time_error_relative_tolerance = 1e-7
  temperature_time_absolute_tolerance = 1e-5
  time_error_safety_factor = 0.9
  restart = study_adaptive_first.checkpoint
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
  exodus = study_adaptive_resume_results.e
  history = study_adaptive_resume_history.csv
  checkpoint = study_adaptive_resume.checkpoint
[]
