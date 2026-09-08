[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]

[Mesh]
  type = exodus
  file = ../moose/m23_pcmi_coupled_cladding_rz_mesh.e
[]

[TimeFunctions]
  [power]
    type = piecewise_linear
    times = 0 20
    values = 0 1
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

    [creep]
      function = norton
      coefficient = 1e-5
      reference_stress = 5e6
      stress_exponent = 3
    []

    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 4e6
      hardening_modulus = 2e9
    []
  []

[]

[Regions]
  [fuel]
    block = fuel
    element = cax4t
    material = fuel
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 2e8
    heat_source_function = power
    heat_source_time_evaluation = interval_average
  []

  [cladding]
    block = clad
    element = cax4t
    material = cladding
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]

[Contact]
  [fuel_cladding]
    primary = clad_left
    secondary = fuel_right

    [thermal]
      gap_conductivity = 0.4
      discretization = node_to_surface
      minimum_gap = 1e-6
    []

    [mechanical]
      discretization = node_to_surface
      sliding = finite
      formulation = penalty
      penalty = 1e14
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

  [fuel_bottom]
    type = dirichlet
    boundary = fuel_bottom
    field = axial_displacement
    value = 0
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
[]

[Executioner]
  type = transient
  end_time = 20
  initial_time_step = 0.0078125
  minimum_time_step = 0.0078125
  maximum_time_step = 0.0078125
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 20
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
  csv = transient_b15_time_00078125_summary.csv
  exodus = transient_b15_time_00078125_results.e
[]
