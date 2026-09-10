[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = ../meshes/b526_friction_reversal.e
[]
[Materials]
  [primary]
    [thermal]
      function = linear_temperature_thermophysical
      conductivity = 15
      density = 100
      specific_heat = 1
      reference_temperature = 300
      conductivity_temperature_coefficient = 0.015
      density_temperature_coefficient = 0
      specific_heat_temperature_coefficient = 0.001
    []
    [elasticity]
      function = linear_temperature_isotropic
      young_modulus = 1.2e8
      poisson_ratio = 0.28
      reference_temperature = 300
      young_modulus_temperature_coefficient = -1e5
      poisson_ratio_temperature_coefficient = 0
    []
    [eigenstrains]
      [thermal]
        function = isotropic_thermal_expansion
        thermal_expansion = 8e-6
        reference_temperature = 300
      []
    []
  []
  [secondary]
    [thermal]
      function = linear_temperature_thermophysical
      conductivity = 10
      density = 100
      specific_heat = 1
      reference_temperature = 300
      conductivity_temperature_coefficient = 0.01
      density_temperature_coefficient = 0
      specific_heat_temperature_coefficient = 0.001
    []
    [elasticity]
      function = linear_temperature_isotropic
      young_modulus = 1e8
      poisson_ratio = 0.3
      reference_temperature = 300
      young_modulus_temperature_coefficient = -5e4
      poisson_ratio_temperature_coefficient = 0
    []
    [eigenstrains]
      [thermal]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 300
      []
    []
  []
[]
[Regions]
  [primary]
    element = c3d8t
    block = primary
    material = primary
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary]
    element = c3d8t
    block = secondary
    material = secondary
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]
[TimeFunctions]
  [secondary_temperature]
    type = piecewise_linear
    times = 0 0.1 0.2 0.3 0.4
    values = 300 400 400 400 400
  []
  [normal_x]
    type = piecewise_linear
    times = 0 0.1 0.2 0.3 0.4
    values = 0 -0.001 -0.001 -0.001 -0.001
  []
  [tangential_y]
    type = piecewise_linear
    times = 0 0.1 0.2 0.3 0.4
    values = 0 0.00002 0.012 -0.004 -0.00398
  []
[]
[Contact]
  [coupled_contact]
    primary = primary_contact
    secondary = secondary_contact
    [thermal]
      law = affine
      conductance = 50
      pressure_derivative = 0.001
    []
    [mechanical]
      formulation = penalty
      penalty = 1e9
      mu = 0.05
      slip_tolerance = 0.005
      discretization = surface_to_surface
      sliding = finite
    []
  []
[]
[BoundaryConditions]
  [primary_temperature]
    type = dirichlet
    boundary = primary_outer
    field = temperature
    value = 300
  []
  [primary_fix_x]
    type = dirichlet
    boundary = primary_outer
    field = displacement_x
    value = 0
  []
  [primary_fix_y]
    type = dirichlet
    boundary = primary_outer
    field = displacement_y
    value = 0
  []
  [primary_fix_z]
    type = dirichlet
    boundary = primary_outer
    field = displacement_z
    value = 0
  []
  [secondary_temperature]
    type = dirichlet
    boundary = secondary_outer
    field = temperature
    value = 1
    function = secondary_temperature
  []
  [secondary_normal_x]
    type = dirichlet
    boundary = secondary_outer
    field = displacement_x
    value = 1
    function = normal_x
  []
  [secondary_tangential_y]
    type = dirichlet
    boundary = secondary_outer
    field = displacement_y
    value = 1
    function = tangential_y
  []
  [secondary_tangential_z]
    type = dirichlet
    boundary = secondary_outer
    field = displacement_z
    value = 0
  []
[]
[Executioner]
  type = transient
  end_time = 0.4
  initial_time_step = 0.02
  minimum_time_step = 0.02
  maximum_time_step = 0.02
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-12
  maximum_iterations = 60
  linear_solver = direct
  preconditioner = lu
  field_residual_scaling = true
  residual_reduction_tolerance = 1e-10
  temperature_residual_absolute_tolerance = 1e-6
  mechanical_residual_absolute_tolerance = 1e-4
[]
[Outputs]
  console = false
  csv = transient_b526_friction_reversal_summary.csv
  exodus = transient_b526_friction_reversal_results.e
[]
