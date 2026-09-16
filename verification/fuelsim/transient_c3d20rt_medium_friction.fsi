[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../moose/m58_integrated_hex20_mesh.e
[]

[TimeFunctions]
  [power]
    type = piecewise_linear
    times = 0 0.2 0.5 0.8 1
    values = 0 0.5 1 1.2 1
  []
  [internal_pressure]
    type = piecewise_linear
    times = 0 0.2 0.5 0.8 1
    values = 0 0.3 1 0.8 1
  []
  [external_pressure]
    type = piecewise_linear
    times = 0 0.2 0.5 0.8 1
    values = 0 0.8 1 1.2 1
  []
  [axial_slide]
    type = piecewise_linear
    times = 0 0.2 0.5 0.8 1
    values = 0 6e-5 3e-4 4.8e-4 6e-4
  []
  [cladding_axial_slide]
    type = piecewise_linear
    times = 0 0.2 0.5 0.8 1
    values = 0 2e-6 1e-5 1.6e-5 2e-5
  []
[]

[Materials]
  [fuel]
    [thermal]
      function = inverse_temperature_thermophysical
      conductivity_inverse_temperature = 3824
      conductivity_constant = 0.61
      density = 10970
      specific_heat = 300
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.316
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 600
      []
    []
  []
  [cladding]
    [thermal]
      function = constant_thermophysical
      conductivity = 16
      density = 6500
      specific_heat = 330
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 7.5e10
      poisson_ratio = 0.3
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 5e-6
        reference_temperature = 600
      []
    []
    [creep]
      function = norton
      coefficient = 1e-10
      reference_stress = 5e5
      stress_exponent = 3
    []
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 1e6
      hardening_modulus = 2e10
    []
  []
[]

[Regions]
  [fuel]
    element = c3d20rt
    block = fuel
    material = fuel
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 2e8
    heat_source_function = power
    heat_source_time_evaluation = interval_average
  []
  [cladding]
    element = c3d20rt
    block = clad
    material = cladding
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]

[Contact]
  [fuel_cladding]
    primary = clad_rmin
    secondary = fuel_outer
    [thermal]
      gap_conductivity = 0.004
      minimum_gap = 1e-6
    []
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      sliding = finite
      penalty = 1e10
      mu = 0.002
    []
  []
[]

[BoundaryConditions]
  [fuel_bottom_x]
    type = dirichlet
    boundary = fuel_bottom
    field = displacement_x
    value = 0
  []
  [fuel_bottom_y]
    type = dirichlet
    boundary = fuel_bottom
    field = displacement_y
    value = 0
  []
  [fuel_top_slide]
    type = dirichlet
    boundary = fuel_top
    field = displacement_z
    value = 1
    function = axial_slide
  []
  [cladding_bottom_x]
    type = dirichlet
    boundary = clad_bottom
    field = displacement_x
    value = 0
  []
  [cladding_bottom_y]
    type = dirichlet
    boundary = clad_bottom
    field = displacement_y
    value = 0
  []
  [cladding_bottom_z]
    type = dirichlet
    boundary = clad_bottom
    field = displacement_z
    value = 0
  []
  [cladding_top_slide]
    type = dirichlet
    boundary = clad_top
    field = displacement_z
    value = 1
    function = cladding_axial_slide
  []
  [cladding_outer_temperature]
    type = dirichlet
    boundary = clad_rmax
    field = temperature
    value = 600
  []
  [cladding_internal_pressure]
    type = pressure
    boundary = clad_rmin
    value = 5e5
    function = internal_pressure
    configuration = current
  []
  [cladding_external_pressure]
    type = pressure
    boundary = clad_rmax
    value = 2e6
    function = external_pressure
    configuration = current
  []
[]

[Executioner]
  type = transient
  end_time = 1
  initial_time_step = 0.05
  minimum_time_step = 0.00625
  maximum_time_step = 0.05
  growth_factor = 2
  cutback_factor = 0.5
  maximum_cutbacks = 3
  load_ramp_time = 0
  use_linear_time_predictor = true
[]

[Solver]
  linear_solver = direct
  preconditioner = lu
  direct_factorization = mumps
  jacobian_lag = 3
  line_search = backtracking
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 40
  backtracking_fallback = true
  field_residual_scaling = true
  temperature_residual_absolute_tolerance = 1e-8
  mechanical_residual_absolute_tolerance = 1e-6
[]

[Outputs]
  console = true
  progress_interval = 1
  exodus = transient_c3d20rt_medium_friction_results.e
  exodus_interval = 1
  csv = transient_c3d20rt_medium_friction_summary.csv
[]
