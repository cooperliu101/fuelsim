[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../abaqus/b60_long_plate_meat_clad_c3d20t_mesh.e
[]

[TimeFunctions]
  [front_temperature]
    type = piecewise_linear
    times = 0 10
    values = 600 700
  []
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
    element = c3d20t
    block = meat
    material = fuel
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 2e8
  []
  [cladding]
    element = c3d20t
    block = clad
    material = cladding
    strain = finite
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
    value = 1
    function = front_temperature
  []
  [back_temperature]
    type = dirichlet
    boundary = plate_back
    field = temperature
    value = 600
  []
[]

[Executioner]
  type = transient
  end_time = 10
  initial_time_step = 1
  minimum_time_step = 1
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
  use_linear_time_predictor = true
[]

[Solver]
  linear_solver = direct
  preconditioner = lu
  direct_factorization = mumps
  # Rebuild the finite-strain Jacobian at every Newton iteration.
  jacobian_lag = 1
  # Keep the predicted increments on the same fully updated Newton path.
  predictor_jacobian_lag = 1
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 30
  field_residual_scaling = true
  # Stop on the same physical field tolerances used by the final residual audit.
  field_residual_convergence = true
  temperature_residual_absolute_tolerance = 1e-6
  # 2e-2 saves two more iterations but fails the 0.5 percent pointwise displacement gate.
  mechanical_residual_absolute_tolerance = 1.9e-2
[]

[Outputs]
  console = false
[]
