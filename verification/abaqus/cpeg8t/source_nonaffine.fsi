[Case]
  version = 3
  problem = transient
  geometry = generalized_plane_strain
[]
[Mesh]
  type = exodus
  file = boundary.e
[]
[TimeFunctions]
  [ramp]
    type = piecewise_linear
    times = 0 1
    values = 0 1
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
  [solid]
    element = cpeg8t
    block = solid
    material = solid
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 1e5
  []
[]
[GeneralizedPlaneStrain]
  [solid]
    blocks = solid
    initial_thickness = 0.1
    u3 = 0
    rotation_x = 0
    rotation_y = 0
  []
[]
[BoundaryConditions]
  [x1]
    type = dirichlet
    boundary = node1
    field = displacement_x
    value = 0
  []
  [x2]
    type = dirichlet
    boundary = node2
    field = displacement_x
    value = 0
  []
  [x3]
    type = dirichlet
    boundary = node3
    field = displacement_x
    value = 0
  []
  [x4]
    type = dirichlet
    boundary = node4
    field = displacement_x
    value = 0
  []
  [x5]
    type = dirichlet
    boundary = node5
    field = displacement_x
    value = 0
  []
  [x6]
    type = dirichlet
    boundary = node6
    field = displacement_x
    value = 0.002
    function = ramp
  []
  [x7]
    type = dirichlet
    boundary = node7
    field = displacement_x
    value = 0
  []
  [x8]
    type = dirichlet
    boundary = node8
    field = displacement_x
    value = -0.001
    function = ramp
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
  initial_time_step = 1
  minimum_time_step = 1
  maximum_time_step = 1
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
  exodus = source_nonaffine_results.e
  history = source_nonaffine_history.csv
[]
