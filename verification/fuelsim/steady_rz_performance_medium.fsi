# Steady CAX4T performance case: explicit production input, SI units.
[Case]
  version = 3
  problem = steady
  geometry = axisymmetric_rz
[]

[Mesh]
  type = exodus
  file = ../meshes/rz_performance_medium.e
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
    element = cax4t
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 2e8
  []
  [cladding]
    block = clad
    material = cladding
    element = cax4t
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
      discretization = node_to_surface
      gap_conductivity = 0.4
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
  type = steady
  load_steps = 20
  maximum_cutbacks = 0
[]

[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-11
  step_tolerance = 1e-16
  maximum_iterations = 100
  line_search = backtracking
  backtracking_fallback = false
  field_residual_scaling = false
  temperature_residual_scale = 1
  mechanical_residual_scale = 100
  temperature_residual_absolute_tolerance = 1e-9
  mechanical_residual_absolute_tolerance = 1e-7
[]

[Outputs]
  console = true
  exodus = steady_rz_performance_medium_results.e
  csv = steady_rz_performance_medium_summary.csv
[]
