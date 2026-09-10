[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../abaqus/b549_small_c3d20t_mesh.e
[]

[TimeFunctions]
  [secondary_temperature]
    type = piecewise_linear
    times = 0 0.4
    values = 300 400
  []
  [outer_pressure]
    type = piecewise_linear
    times = 0 0.4
    values = 0 175000
  []
  [tangential_y]
    type = piecewise_linear
    times = 0 0.4
    values = 0 0
  []
  [tangential_z]
    type = piecewise_linear
    times = 0 0.4
    values = 0 0
  []
[]

[Materials]
  [primary]
    [thermal]
      function = constant_thermophysical
      conductivity = 15
      density = 100
      specific_heat = 1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1.2e8
      poisson_ratio = 0.28
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 8e-6
        reference_temperature = 300
      []
    []
  []
  [secondary]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 100
      specific_heat = 1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e8
      poisson_ratio = 0.3
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 300
      []
    []
    [creep]
      function = norton
      coefficient = 1e-4
      reference_stress = 2e5
      stress_exponent = 3
    []
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 1e5
      hardening_modulus = 1e7
    []
  []
[]

[Regions]
  [primary]
    element = c3d20t
    block = primary
    material = primary
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary]
    element = c3d20t
    block = secondary
    material = secondary
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
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
  [primary_contact_temperature]
    type = dirichlet
    boundary = primary_contact
    field = temperature
    value = 1
    function = secondary_temperature
  []
  [secondary_contact_temperature]
    type = dirichlet
    boundary = secondary_contact
    field = temperature
    value = 300
  []
  [secondary_contact_fix_x]
    type = dirichlet
    boundary = secondary_contact
    field = displacement_x
    value = 0
  []
  [secondary_temperature]
    type = dirichlet
    boundary = secondary_outer
    field = temperature
    value = 1
    function = secondary_temperature
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
    value = 1
    function = tangential_z
  []
  [secondary_outer_pressure]
    type = pressure
    boundary = secondary_outer
    value = 1
    function = outer_pressure
    configuration = current
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
  linear_solver = direct
  preconditioner = lu
  direct_factorization = mumps
  line_search = backtracking
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-12
  step_tolerance = 1e-12
  maximum_iterations = 60
  field_residual_scaling = true
  residual_reduction_tolerance = 1e-10
  temperature_residual_absolute_tolerance = 1e-6
  mechanical_residual_absolute_tolerance = 1e-4
[]

[Outputs]
  console = false
  csv = transient_b549_small_c3d20t_summary.csv
  exodus = transient_b549_small_c3d20t_results.e
  exodus_interval = 20
  checkpoint = transient_b549_small_c3d20t.checkpoint
  checkpoint_interval = 20
[]
