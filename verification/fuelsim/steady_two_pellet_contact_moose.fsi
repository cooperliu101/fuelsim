[Case]
  version = 2
  problem = steady
[]

[Mesh]
  type = exodus
  file = ../moose/m33_two_pellet_contact_rz_mesh.e
[]

[Materials]
  [lower]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1
      specific_heat = 1
    []

    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.3
    []

    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 600
      []
    []
  []

  [upper]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1
      specific_heat = 1
    []

    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.3
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
  [lower]
    block = lower_pellet
    material = lower
    strain = small
    initial_temperature = 800
    volumetric_heat_source = 0
  []
  [upper]
    block = upper_pellet
    material = upper
    strain = small
    initial_temperature = 800
    volumetric_heat_source = 0
  []
[]

[Contact]
  [pellet_stack]
    primary = upper_bottom
    secondary = lower_top
    [thermal]
      gap_conductivity = 0.2
      minimum_gap = 1e-6
    []
    [mechanical]
      formulation = penalty
      penalty = 1e14
    []
  []
[]

[BoundaryConditions]
  [lower_axis]
    type = dirichlet
    boundary = lower_left
    field = radial_displacement
    value = 0
  []
  [upper_axis]
    type = dirichlet
    boundary = upper_left
    field = radial_displacement
    value = 0
  []
  [lower_bottom]
    type = dirichlet
    boundary = lower_bottom
    field = axial_displacement
    value = 0
  []
  [upper_top]
    type = dirichlet
    boundary = upper_top
    field = axial_displacement
    value = 0
  []
  [lower_temperature]
    type = dirichlet
    boundary = lower_right
    field = temperature
    value = 800
  []
  [upper_temperature]
    type = dirichlet
    boundary = upper_right
    field = temperature
    value = 800
  []
[]

[Executioner]
  type = steady
  load_steps = 1
[]

[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 50
[]

[Outputs]
  console = false
[]
