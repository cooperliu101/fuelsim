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
  [extension]
    type = piecewise_linear
    times = 0 1 2
    values = 0 1 1
  []
  [x1]
    type = piecewise_linear
    times = 0 1 2
    values = 0 -0.00025 0.01495
  []
  [x2]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.00025 -0.00505
  []
  [x3]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.00015 -0.01505
  []
  [x4]
    type = piecewise_linear
    times = 0 1 2
    values = 0 -0.00015 0.00495
  []
  [x5]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0 0.005
  []
  [x6]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.0002 -0.01005
  []
  [x7]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0 -0.005
  []
  [x8]
    type = piecewise_linear
    times = 0 1 2
    values = 0 -0.0002 0.00995
  []
  [y1]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.00005 -0.00525
  []
  [y2]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.00005 0.01525
  []
  [y3]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.00005 0.00515
  []
  [y4]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.00005 -0.01515
  []
  [y5]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0 0.005
  []
  [y6]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.00005 0.0102
  []
  [y7]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0 -0.005
  []
  [y8]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.00005 -0.0102
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
    volumetric_heat_source = 0
  []
[]
[GeneralizedPlaneStrain]
  [solid]
    blocks = solid
    initial_thickness = 0.1
    u3 = 0.001
    u3_function = extension
    rotation_x = 0
    rotation_y = 0
  []
[]
[BoundaryConditions]
  [x1]
    type = dirichlet
    boundary = node1
    field = displacement_x
    value = 1
    function = x1
  []
  [x2]
    type = dirichlet
    boundary = node2
    field = displacement_x
    value = 1
    function = x2
  []
  [x3]
    type = dirichlet
    boundary = node3
    field = displacement_x
    value = 1
    function = x3
  []
  [x4]
    type = dirichlet
    boundary = node4
    field = displacement_x
    value = 1
    function = x4
  []
  [x5]
    type = dirichlet
    boundary = node5
    field = displacement_x
    value = 1
    function = x5
  []
  [x6]
    type = dirichlet
    boundary = node6
    field = displacement_x
    value = 1
    function = x6
  []
  [x7]
    type = dirichlet
    boundary = node7
    field = displacement_x
    value = 1
    function = x7
  []
  [x8]
    type = dirichlet
    boundary = node8
    field = displacement_x
    value = 1
    function = x8
  []
  [y1]
    type = dirichlet
    boundary = node1
    field = displacement_y
    value = 1
    function = y1
  []
  [y2]
    type = dirichlet
    boundary = node2
    field = displacement_y
    value = 1
    function = y2
  []
  [y3]
    type = dirichlet
    boundary = node3
    field = displacement_y
    value = 1
    function = y3
  []
  [y4]
    type = dirichlet
    boundary = node4
    field = displacement_y
    value = 1
    function = y4
  []
  [y5]
    type = dirichlet
    boundary = node5
    field = displacement_y
    value = 1
    function = y5
  []
  [y6]
    type = dirichlet
    boundary = node6
    field = displacement_y
    value = 1
    function = y6
  []
  [y7]
    type = dirichlet
    boundary = node7
    field = displacement_y
    value = 1
    function = y7
  []
  [y8]
    type = dirichlet
    boundary = node8
    field = displacement_y
    value = 1
    function = y8
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
  end_time = 2
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
  exodus = bending_rotation_results.e
  history = bending_rotation_history.csv
[]
