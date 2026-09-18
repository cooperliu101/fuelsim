[Case]
  version = 3
  physics = thermal
  problem = steady
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = dcax8.e
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
    element = dcax8
    material = solid
    initial_temperature = 300
    volumetric_heat_source = 100000
  []
[]
[BoundaryConditions]
  [hot]
    type = heat_flux
    boundary = left
    value = 1000
  []
  [cold]
    type = convection
    boundary = right
    heat_transfer_coefficient = 100
    ambient_temperature = 300
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
  exodus = dcax8_steady_results.e
[]
