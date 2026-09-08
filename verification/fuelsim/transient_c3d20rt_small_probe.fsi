[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = ../meshes/c3d20rt_probe.e
[]
[TimeFunctions]
  [t1]
    type = piecewise_linear
    times = 0 1 3
    values = 300 301 301
  []
  [t2]
    type = piecewise_linear
    times = 0 1 3
    values = 300 302 302
  []
  [t3]
    type = piecewise_linear
    times = 0 1 3
    values = 300 303 303
  []
  [t4]
    type = piecewise_linear
    times = 0 1 3
    values = 300 304 304
  []
  [t5]
    type = piecewise_linear
    times = 0 1 3
    values = 300 305 305
  []
  [t6]
    type = piecewise_linear
    times = 0 1 3
    values = 300 306 306
  []
  [t7]
    type = piecewise_linear
    times = 0 1 3
    values = 300 307 307
  []
  [t8]
    type = piecewise_linear
    times = 0 1 3
    values = 300 308 308
  []
  [power]
    type = piecewise_linear
    times = 0 2 3
    values = 0 0 1
  []
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 2
      specific_heat = 3
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 200000
      poisson_ratio = 0.25
    []
  []
[]
[Regions]
  [solid]
    block = solid
    material = solid
    element = c3d20rt
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 10
    heat_source_function = power
  []
[]
[BoundaryConditions]
  [uy]
    type = dirichlet
    boundary = all
    field = displacement_y
    value = 0
  []
  [uz]
    type = dirichlet
    boundary = all
    field = displacement_z
    value = 0
  []
  [uxzero]
    type = dirichlet
    boundary = xzero
    field = displacement_x
    value = 0
  []
  [uxhalf]
    type = dirichlet
    boundary = xhalf
    field = displacement_x
    value = 0.15
    scale_with_load = true
  []
  [uxone]
    type = dirichlet
    boundary = xone
    field = displacement_x
    value = 0.2
    scale_with_load = true
  []
  [t1]
    type = dirichlet
    boundary = n1
    field = temperature
    value = 1
    function = t1
  []
  [t2]
    type = dirichlet
    boundary = n2
    field = temperature
    value = 1
    function = t2
  []
  [t3]
    type = dirichlet
    boundary = n3
    field = temperature
    value = 1
    function = t3
  []
  [t4]
    type = dirichlet
    boundary = n4
    field = temperature
    value = 1
    function = t4
  []
  [t5]
    type = dirichlet
    boundary = n5
    field = temperature
    value = 1
    function = t5
  []
  [t6]
    type = dirichlet
    boundary = n6
    field = temperature
    value = 1
    function = t6
  []
  [t7]
    type = dirichlet
    boundary = n7
    field = temperature
    value = 1
    function = t7
  []
  [t8]
    type = dirichlet
    boundary = n8
    field = temperature
    value = 1
    function = t8
  []
[]
[Executioner]
  type = transient
  end_time = 3
  initial_time_step = 1
  minimum_time_step = 1
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 1
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  maximum_iterations = 20
[]
[Outputs]
  console = true
  exodus = transient_c3d20rt_small_probe_results.e
[]
