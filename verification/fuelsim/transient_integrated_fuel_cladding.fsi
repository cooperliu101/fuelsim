[Case]
  version = 1
  problem = transient
[]

[Mesh]
  type = exodus
  file = ../moose/m57_integrated_fuel_cladding_rz_mesh.e
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

[Regions]
  [fuel]
    block = fuel
    strain = finite
    conductivity_inverse_temperature = 3824
    conductivity_constant = 0.61
    young_modulus = 2e11
    poisson_ratio = 0.316
    thermal_expansion = 1e-5
    reference_temperature = 600
    initial_temperature = 600
    volumetric_heat_source = 2e8
    heat_source_function = power
    density = 10970
    specific_heat = 300
    inelastic_model = elastic
  []

  [cladding]
    block = clad
    strain = finite
    conductivity_inverse_temperature = 0
    conductivity_constant = 16
    young_modulus = 7.5e10
    poisson_ratio = 0.3
    thermal_expansion = 5e-6
    reference_temperature = 600
    initial_temperature = 600
    volumetric_heat_source = 0
    density = 6500
    specific_heat = 330
    inelastic_model = norton_creep_j2_plasticity
    creep_coefficient = 1e-10
    creep_reference_stress = 5e5
    creep_exponent = 3
    yield_stress = 1e6
    hardening_modulus = 2e10
  []
[]

[Contact]
  [fuel_cladding]
    primary = clad_left
    secondary = fuel_right

    [thermal]
      gap_conductivity = 0.004
      minimum_gap = 1e-6
    []

    [mechanical]
      formulation = penalty
      penalty = 1e12
      mu = 0.002
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
  []
  [cladding_external_pressure]
    type = pressure
    boundary = clad_right
    value = 2e6
    function = external_pressure
  []
[]

[Executioner]
  type = transient
  end_time = 6
  initial_time_step = 1.5
  minimum_time_step = 0.03125
  maximum_time_step = 1.5
  growth_factor = 2
  cutback_factor = 0.5
  maximum_cutbacks = 8
  load_ramp_time = 0
  target_nonlinear_iterations = 8
  iteration_window = 2
  time_error_relative_tolerance = 1e-2
  temperature_time_absolute_tolerance = 1e-2
  displacement_time_absolute_tolerance = 1e-7
  strain_history_time_absolute_tolerance = 1e-4
  stress_history_time_absolute_tolerance = 1e7
  time_error_safety_factor = 0.9
[]

[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 100
  backtracking_fallback = true
  field_residual_scaling = true
  temperature_residual_absolute_tolerance = 1e-8
  mechanical_residual_absolute_tolerance = 1e-4
[]

[Outputs]
  console = false
[]
