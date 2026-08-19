[Case]
  version = 3
  problem = steady
  geometry = axisymmetric_rz
[]

[Mesh]
  type = exodus
  file = ../moose/m0_simple_fuel_rz_mesh.e
[]

[Materials]
  [fuel]
    [thermal]
      function = inverse_temperature_thermophysical
      conductivity_inverse_temperature = 3824
      conductivity_constant = 0.61
      density = 1
      specific_heat = 1
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

[]

[Regions]
  [fuel]
    block_id = 0
    material = fuel
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 2e8
  []
[]

[BoundaryConditions]
  [axis]
    type = dirichlet
    boundary = fuel_left
    field = radial_displacement
    value = 0
  []
  [bottom]
    type = dirichlet
    boundary = fuel_bottom
    field = axial_displacement
    value = 0
  []
  [outer_temperature]
    type = dirichlet
    boundary = fuel_right
    field = temperature
    value = 600
  []
[]

[Executioner]
  type = steady
  load_steps = 1
[]

[Solver]
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 40
[]

[Outputs]
  console = false
[]
