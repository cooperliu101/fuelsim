# Finite-strain CAX8RT: 20 equilibrium increments with committed material history.
# Thermal storage is disabled, so every increment solves steady heat balance.
[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]

[Mesh]
  type = exodus
  file = ../meshes/rz_performance_medium_cax8t.e
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
  []
[]

[Regions]
  [fuel]
    block = fuel
    material = fuel
    element = cax8rt
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 2e8
  []
  [cladding]
    block = clad
    material = cladding
    element = cax8rt
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
      discretization = node_to_surface
      gap_conductivity = 0.4
      minimum_gap = 1e-6
    []
    [mechanical]
      discretization = node_to_surface
      sliding = finite
      formulation = penalty
      penalty = 1e14
      mu = 0.2
      slip_tolerance = 0.001
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
  initial_time_step = 1
  minimum_time_step = 1
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 20
  include_thermal_time_term = false
  use_linear_time_predictor = true
[]

[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-11
  step_tolerance = 1e-16
  maximum_iterations = 100
  line_search = basic
  backtracking_fallback = true
  field_residual_scaling = false
  temperature_residual_scale = 1
  mechanical_residual_scale = 100
  temperature_residual_absolute_tolerance = 1e-9
  mechanical_residual_absolute_tolerance = 1e-7
[]

[Outputs]
  console = true
[]
