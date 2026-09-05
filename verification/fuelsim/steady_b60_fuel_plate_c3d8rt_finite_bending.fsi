[Case]
  version = 3
  problem = steady
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../moose/b60_long_plate_meat_clad_mesh.e
[]

[Materials]
  [fuel]
    [thermal]
      function = constant_thermophysical
      conductivity = 3
      density = 10970
      specific_heat = 300
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.30
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-4
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
      young_modulus = 1e11
      poisson_ratio = 0.32
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
    block = meat
    material = fuel
    strain = finite
    element = c3d8rt
    initial_temperature = 600
    volumetric_heat_source = 2e8
  []
  [cladding]
    block = clad
    material = cladding
    strain = finite
    element = c3d8rt
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]

[BoundaryConditions]
  [clamp]
    type = dirichlet
    boundary = plate_left
    field = displacement_x
    value = 0
  []
  [clamp_y]
    type = dirichlet
    boundary = plate_left
    field = displacement_y
    value = 0
  []
  [clamp_z]
    type = dirichlet
    boundary = plate_left
    field = displacement_z
    value = 0
  []
  [front_temperature]
    type = dirichlet
    boundary = plate_front
    field = temperature
    value = 700
  []
  [back_temperature]
    type = dirichlet
    boundary = plate_back
    field = temperature
    value = 600
  []
[]

[Executioner]
  type = steady
  load_steps = 1
  # Use the full-load small-strain equilibrium only as the finite-strain initial guess.
  use_small_strain_predictor = true
[]

[Solver]
  linear_solver = direct
  preconditioner = lu
  direct_factorization = mumps
  jacobian_lag = 1
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 30
  field_residual_scaling = true
  # Stop on the same physical field tolerances used by the final residual audit.
  field_residual_convergence = true
  temperature_residual_absolute_tolerance = 1e-6
  mechanical_residual_absolute_tolerance = 1e-3
[]

[Outputs]
  console = false
  csv = steady_b60_fuel_plate_c3d8rt_finite_bending_summary.csv
  exodus = steady_b60_fuel_plate_c3d8rt_finite_bending_results.e
[]
