[Case]
  version = 3
  problem = transient
  geometry = generalized_plane_strain
[]
[Mesh]
  type = exodus
  file = contact.e
[]
[TimeFunctions]
  [ramp]
    type = piecewise_linear
    times = 0 1
    values = 0 1
  []
  [delayed]
    type = piecewise_linear
    times = 0 0.5 1
    values = 0 0 1
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
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [upper]
    element = cpeg8t
    block = upper
    material = solid
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]
[GeneralizedPlaneStrain]
  [lower]
    blocks = lower
    initial_thickness = 0.1
    u3 = 0.001
    u3_function = ramp
    rotation_x = 0.02
    rotation_x_function = ramp
    rotation_y = -0.03
    rotation_y_function = ramp
  []
  [upper]
    blocks = upper
    initial_thickness = 0.2
    rotation_x = -0.04
    rotation_x_function = delayed
    rotation_y = 0.05
    rotation_y_function = delayed
  []
[]
[BoundaryConditions]
  [fix_x]
    type = dirichlet
    boundary = field
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = field
    field = displacement_y
    value = 0
  []
  [temperature]
    type = dirichlet
    boundary = field
    field = temperature
    value = 300
  []
[]
[Executioner]
  type = transient
  end_time = 1
  initial_time_step = 0.5
  minimum_time_step = 0.5
  maximum_time_step = 0.5
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-11
  step_tolerance = 1e-14
  maximum_iterations = 20
  linear_solver = direct
  direct_factorization = mumps
[]
[Outputs]
  console = true
  exodus = sections_results.e
  history = sections_history.csv
  checkpoint = sections.chk
[]
