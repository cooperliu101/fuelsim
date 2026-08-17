[Case]
  version = 3
  problem = steady
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../moose/b35_hex8_shared_meat_clad_mesh.e
[]

[Materials]
  [meat]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 10970
      specific_heat = 300
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
  [cladding]
    [thermal]
      function = constant_thermophysical
      conductivity = 20
      density = 6500
      specific_heat = 330
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 5e9
      poisson_ratio = 0.3
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 5e-6
        reference_temperature = 300
      []
    []
  []
[]

[Regions]
  [meat]
    block = meat
    material = meat
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [clad]
    block = clad
    material = cladding
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]

[Contact]
[]

[BoundaryConditions]
  [temperature_left]
    type = dirichlet
    boundary = plate_left
    field = temperature
    value = 300
  []
  [temperature_right]
    type = dirichlet
    boundary = plate_right
    field = temperature
    value = 600
  []
  [fix_x]
    type = dirichlet
    boundary = plate_left
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = plate_left
    field = displacement_y
    value = 0
  []
  [fix_z]
    type = dirichlet
    boundary = plate_left
    field = displacement_z
    value = 0
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
