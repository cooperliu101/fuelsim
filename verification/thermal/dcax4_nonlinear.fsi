[Case]
  version = 3
  physics = thermal
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = dcax4.e
[]
[Materials]
  [solid]
    [thermal]
      function = linear_temperature_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
      reference_temperature = 300
      conductivity_temperature_coefficient = 0.05
      specific_heat_temperature_coefficient = 1
    []
  []
[]
[Regions]
  [solid]
    block = solid
    element = dcax4
    material = solid
    initial_temperature = 300
    volumetric_heat_source = 0
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
  end_time = 20
  initial_time_step = 2
  minimum_time_step = 2
  maximum_time_step = 2
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
  exodus = dcax4_nonlinear_results.e
  history = dcax4_nonlinear_history.csv
[]
