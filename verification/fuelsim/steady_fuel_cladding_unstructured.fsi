[Case]
  version = 2
  problem = steady
[]

[Mesh]
  type = exodus
  file = ../moose/m1_fuel_cladding_unstructured_rz_mesh.e
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

  [cladding]
    [thermal]
      function = constant_thermophysical
      conductivity = 16
      density = 1
      specific_heat = 1
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
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 2e8
  []

  [cladding]
    block = clad
    material = cladding
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
      gap_conductivity = 0.4
      minimum_gap = 1e-6
    []

    [mechanical]
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
[]

[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 50
[]

[Outputs]
  console = true
[]
