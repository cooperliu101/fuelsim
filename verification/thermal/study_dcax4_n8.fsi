[Case]
  version = 3
  physics = thermal
  problem = steady
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = study_dcax4_n8.e
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
  [hot]
    type = dirichlet
    boundary = left
    field = temperature
    value = 400
  []
  [cold]
    type = dirichlet
    boundary = right
    field = temperature
    value = 300
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
  maximum_iterations = 40
[]
[Outputs]
  console = true
  exodus = study_dcax4_n8_results.e
[]
