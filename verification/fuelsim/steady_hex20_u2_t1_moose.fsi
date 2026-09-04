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
      density = 6000
      specific_heat = 1000
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.25
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 300
      []
    []
  []
[]

[Regions]
  [solid]
    block_id = 0
    material = solid
    strain = small
    initial_temperature = 400
    volumetric_heat_source = 0
  []
[]

[BoundaryConditions]
  [temperature]
    type = dirichlet
    boundary = left
    field = temperature
    value = 400
  []
  [temperature_right]
    type = dirichlet
    boundary = right
    field = temperature
    value = 400
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
  load_steps = 1
[]

[Solver]
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-11
  step_tolerance = 1e-12
  maximum_iterations = 20
  linear_solver = direct
  preconditioner = lu
[]

[Outputs]
  console = false
  exodus = steady_hex20_u2_t1_moose_results.e
[]
