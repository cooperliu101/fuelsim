[Case]
  version = 3
  problem = steady
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = ../moose/b6_hex20_u2_t1.e
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 2000
      specific_heat = 1000
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e6
      poisson_ratio = 0.25
    []
    [eigenstrains]
      [expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-3
        reference_temperature = 300
      []
    []
  []
[]
[Regions]
  [solid]
    block_id = 0
    material = solid
    element = c3d20rt
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
    body_acceleration = 9.81 0 0
  []
[]
[BoundaryConditions]
  [left_temperature]
    type = dirichlet
    boundary = left
    field = temperature
    value = 500
  []
  [right_temperature]
    type = dirichlet
    boundary = right
    field = temperature
    value = 500
  []
  [fix_x]
    type = dirichlet
    boundary = left
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = bottom
    field = displacement_y
    value = 0
  []
  [fix_z]
    type = dirichlet
    boundary = back
    field = displacement_z
    value = 0
  []
[]
[Executioner]
  type = steady
  load_steps = 4
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-11
  step_tolerance = 1e-12
  maximum_iterations = 40
[]
[Outputs]
  console = false
  csv = steady_gravity_hex20_summary.csv
  exodus = steady_gravity_hex20_results.e
[]
