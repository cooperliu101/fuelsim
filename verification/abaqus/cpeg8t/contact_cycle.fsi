[Case]
  version = 3
  problem = transient
  geometry = generalized_plane_strain
[]
[Mesh]
  type = exodus
  file = sliding.e
[]
[TimeFunctions]
  [ramp]
    type = piecewise_linear
    times = 0 3
    values = 0 1
  []
  [vertical]
    type = piecewise_linear
    times = 0 0.5 1 1.5 2 2.5 3
    values = 0 -0.0002 -0.0002 0 0 -0.0002 -0.0002
  []
  [temperature]
    type = piecewise_linear
    times = 0 3
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
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e6
      poisson_ratio = 0.25
    []
  []
[]
[Regions]
  [lower]
    element = cpeg8t
    block = lower
    material = solid
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [upper]
    element = cpeg8t
    block = upper
    material = solid
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]
[GeneralizedPlaneStrain]
  [lower]
    blocks = lower
    initial_thickness = 0.1
    u3 = 0
    rotation_x = 0
    rotation_y = 0
  []
  [upper]
    blocks = upper
    initial_thickness = 0.1
    u3 = 0.01
    u3_function = ramp
    rotation_x = 0.02
    rotation_x_function = ramp
    rotation_y = -0.03
    rotation_y_function = ramp
  []
[]
[Contact]
  [interface]
    primary = primary
    secondary = secondary
    [thermal]
      law = affine
      conductance = 1000
    []
    [mechanical]
      formulation = penalty
      penalty = 1e9
      discretization = surface_to_surface
      sliding = finite
    []
  []
[]
[BoundaryConditions]
  [lower_x]
    type = dirichlet
    boundary = lower
    field = displacement_x
    value = 0
  []
  [lower_y]
    type = dirichlet
    boundary = lower
    field = displacement_y
    value = 0
  []
  [upper_x]
    type = dirichlet
    boundary = upper
    field = displacement_x
    value = 0.015
    function = ramp
  []
  [upper_y]
    type = dirichlet
    boundary = upper
    field = displacement_y
    value = 1
    function = vertical
  []
  [hot]
    type = dirichlet
    boundary = upper
    field = temperature
    value = 1
    function = temperature
  []
  [cold]
    type = dirichlet
    boundary = lower
    field = temperature
    value = 300
  []
[]
[Executioner]
  type = transient
  end_time = 3
  initial_time_step = 0.25
  minimum_time_step = 0.25
  maximum_time_step = 0.25
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
  include_thermal_time_term = true
[]
[Solver]
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-11
  step_tolerance = 1e-14
  maximum_iterations = 30
  linear_solver = direct
  direct_factorization = mumps
[]
[Outputs]
  console = true
  exodus = contact_cycle_results.e
  history = contact_cycle_history.csv
  checkpoint = contact_cycle.chk
[]
