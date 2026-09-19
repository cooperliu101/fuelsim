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
  []
[]
[BoundaryConditions]
[]
[Executioner]
  type = transient
  end_time = 2
  initial_time_step = 0.4
  minimum_time_step = 0.4
  maximum_time_step = 0.4
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
  exodus = study_time_coarse_results.e
  history = study_time_coarse_history.csv
[]
