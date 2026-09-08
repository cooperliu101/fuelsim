[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]

[Mesh]
  type = exodus
  file = ../meshes/b13_integrated_cax8rt.e
[]

[TimeFunctions]
  [power]
    type = piecewise_linear
    times = 0 0.5 1.5 3 4.5 6
    values = 0 0.4 1 0.8 1.2 1
  []
  [internal_pressure]
    type = piecewise_linear
    times = 0 0.5 1.5 3 4.5 6
    values = 0 0.25 1 0.7 1.1 1
  []
  [external_pressure]
    type = piecewise_linear
    times = 0 0.5 1.5 3 4.5 6
    values = 0 0.8 1 1.2 0.9 1
  []
  [axial_slide]
    type = piecewise_linear
    times = 0 0.5 1.5 3 4.5 6
    values = 0 6e-5 2.4e-4 3.6e-4 4.8e-4 6e-4
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
    block = fuel
    element = cax8rt
    material = fuel
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 2e7
    heat_source_function = power
    heat_source_time_evaluation = interval_average
  []

  [cladding]
    block = clad
    element = cax8rt
    material = cladding
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]

[Contact]
  [fuel_cladding]
    primary = clad_left
    secondary = fuel_right

    [thermal]
      gap_conductivity = 0.004
      discretization = node_to_surface
      minimum_gap = 1e-6
    []

    [mechanical]
      discretization = node_to_surface
      sliding = finite
      formulation = penalty
      penalty = 1e12
      mu = 0.002
      slip_tolerance = 0.005
    []
  []
[]

[BoundaryConditions]
  [fuel_axis]
    type = dirichlet
    boundary = fuel_left
    field = radial_displacement
    value = 0
  []
  [fuel_top_slide]
    type = dirichlet
    boundary = fuel_top
    field = axial_displacement
    value = 1
    function = axial_slide
  []
  [cladding_bottom]
    type = dirichlet
    boundary = clad_bottom
    field = axial_displacement
    value = 0
  []
  [cladding_outer_temperature]
    type = dirichlet
    boundary = clad_right
    field = temperature
    value = 600
  []
  [cladding_internal_pressure]
    type = pressure
    boundary = clad_left
    value = 5e5
    function = internal_pressure
    configuration = current
  []
  [cladding_external_pressure]
    type = pressure
    boundary = clad_right
    value = 2e6
    function = external_pressure
    configuration = current
  []
[]

[Executioner]
  type = transient
  end_time = 6
  initial_time_step = .0625
  minimum_time_step = .0625
  maximum_time_step = .0625
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-11
  step_tolerance = 1e-16
  maximum_iterations = 100
  line_search = backtracking
  backtracking_fallback = true
  field_residual_scaling = false
  temperature_residual_scale = 1
  mechanical_residual_scale = 100
  temperature_residual_absolute_tolerance = 1e-9
  mechanical_residual_absolute_tolerance = 1e-7
[]
[Outputs]
  console = true
  csv = transient_b13_integrated_cax8rt_summary.csv
  exodus = transient_b13_integrated_cax8rt_results.e
[]
