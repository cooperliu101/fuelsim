[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = ../meshes/b59_distorted.e
[]
[TimeFunctions]
  [bottom_temperature]
    type = piecewise_linear
    times = 0 10000000 20000000 30000000 40000000
    values = 300 480 600 600 330
  []
[]
[Materials]
  [lower]
    [thermal]
      function = constant_thermophysical
      conductivity = 8
      density = 10000
      specific_heat = 300
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1.2e11
      poisson_ratio = 0.28
    []
    [eigenstrains]
      [thermal]
        function = isotropic_thermal_expansion
        thermal_expansion = 1.8e-5
        reference_temperature = 300
      []
    []
  []
  [upper]
    [thermal]
      function = constant_thermophysical
      conductivity = 24
      density = 6500
      specific_heat = 500
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.3
    []
    [eigenstrains]
      [thermal]
        function = isotropic_thermal_expansion
        thermal_expansion = 7e-6
        reference_temperature = 300
      []
    []
  []
[]
[Regions]
  [lower]
    block = lower
    material = lower
    strain = small
    element = c3d8t
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [upper]
    block = upper
    material = upper
    strain = small
    element = c3d8t
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]
[BoundaryConditions]
  [bottom_temperature]
    type = dirichlet
    boundary = bottom
    field = temperature
    value = 1
    function = bottom_temperature
  []
  [top_temperature]
    type = dirichlet
    boundary = top
    field = temperature
    value = 300
  []
  [fix_x]
    type = dirichlet
    boundary = left_lower
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = left_lower
    field = displacement_y
    value = 0
  []
  [fix_z]
    type = dirichlet
    boundary = left_lower
    field = displacement_z
    value = 0
  []
[]
[Executioner]
  type = transient
  end_time = 40000000
  initial_time_step = 10000000
  minimum_time_step = 10000000
  maximum_time_step = 10000000
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-13
  maximum_iterations = 20
  linear_solver = direct
  preconditioner = lu
  field_residual_scaling = true
  residual_reduction_tolerance = 1e-11
  temperature_residual_absolute_tolerance = 1e-8
  mechanical_residual_absolute_tolerance = 1e-2
[]
[Outputs]
  console = false
  csv = transient_b59_distorted_summary.csv
  exodus = transient_b59_distorted_results.e
[]
