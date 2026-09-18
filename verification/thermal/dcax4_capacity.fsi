[Case]
  version = 3
  physics = thermal
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = dcax4_capacity.e
[]
[TimeFunctions]
  [ramp]
    type = piecewise_linear
    times = 0 1
    values = 300 400
  []
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
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
  [uniform]
    type = dirichlet
    boundary = all
    field = temperature
    value = 1
    function = ramp
  []
[]
[Executioner]
  type = transient
  end_time = 1
  initial_time_step = 1
  minimum_time_step = 1
  maximum_time_step = 1
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
  exodus = dcax4_capacity_results.e
  history = dcax4_capacity_history.csv
[]
