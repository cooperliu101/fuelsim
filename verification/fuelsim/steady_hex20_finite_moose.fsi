[Case]
  version = 3
  problem = steady
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = ../moose/h20_10_finite.e
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
    element = c3d20t
    block_id = 0
    material = solid
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]
[BoundaryConditions]
  [temperature]
    type = dirichlet
    boundary = left
    field = temperature
    value = 300
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
  [traction]
    type = traction
    boundary = right
    field = displacement_x
    value = 1000000
    configuration = current
  []
[]
[Executioner]
  type = steady
  load_steps = 1
[]
[Solver]
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-12
  step_tolerance = 1e-12
  maximum_iterations = 30
[]
[Outputs]
  console = false
  exodus = steady_hex20_finite_moose_results.e
[]
