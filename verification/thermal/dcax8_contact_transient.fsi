[Case]
  version = 3
  physics = thermal
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = dcax8_contact.e
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
  [inner]
    block = inner
    element = dcax8
    material = solid
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [outer]
    block = outer
    element = dcax8
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
  exodus = dcax8_contact_transient_results.e
  history = dcax8_contact_transient_history.csv
  checkpoint = dcax8_contact_transient.checkpoint
[]
