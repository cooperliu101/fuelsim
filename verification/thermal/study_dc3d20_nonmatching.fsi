[Case]
  version = 3
  physics = thermal
  problem = steady
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = study_dc3d20_nonmatching.e
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
  [hot]
    block = hot
    element = dc3d20
    material = solid
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [cold]
    block = cold
    element = dc3d20
    material = solid
    initial_temperature = 300
    volumetric_heat_source = 0
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
  [hot0]
    type = dirichlet
    boundary = hot0
    field = temperature
    value = 400
  []
  [hot1]
    type = dirichlet
    boundary = hot1
    field = temperature
    value = 425
  []
  [hot2]
    type = dirichlet
    boundary = hot2
    field = temperature
    value = 450
  []
  [hot3]
    type = dirichlet
    boundary = hot3
    field = temperature
    value = 475
  []
  [hot4]
    type = dirichlet
    boundary = hot4
    field = temperature
    value = 500
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
  exodus = study_dc3d20_nonmatching_results.e
[]
