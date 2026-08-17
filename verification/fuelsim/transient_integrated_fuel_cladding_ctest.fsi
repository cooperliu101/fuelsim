[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]

[Mesh]
  type = exodus
  file = ../moose/m57_integrated_fuel_cladding_rz_mesh.e
[]

[TimeFunctions]
  [power]
    type = piecewise_linear
    times = 0 0.0833333333333333 0.25 0.5 0.75 1
    values = 0 0.4 1 0.8 1.2 1
  []
  [internal_pressure]
    type = piecewise_linear
    times = 0 0.0833333333333333 0.25 0.5 0.75 1
    values = 0 0.25 1 0.7 1.1 1
  []
  [external_pressure]
    type = piecewise_linear
    times = 0 0.0833333333333333 0.25 0.5 0.75 1
    values = 0 0.8 1 1.2 0.9 1
  []
  [axial_slide]
    type = piecewise_linear
    times = 0 0.0833333333333333 0.25 0.5 0.75 1
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
    material = fuel
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 2e7
    heat_source_function = power
  []

  [cladding]
    block = clad
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
  # The default CTest reference retains the complete load path in one sixth
  # of its physical duration.  The full 6 s path remains the manual MPI
  # benchmark in transient_integrated_fuel_cladding.fsi.
  end_time = 1
  initial_time_step = 0.03125
  minimum_time_step = 0.015625
  maximum_time_step = 0.0625
  growth_factor = 2
  cutback_factor = 0.5
  maximum_cutbacks = 8
  load_ramp_time = 0
  target_nonlinear_iterations = 8
  iteration_window = 2
  time_error_relative_tolerance = 1e-1
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
