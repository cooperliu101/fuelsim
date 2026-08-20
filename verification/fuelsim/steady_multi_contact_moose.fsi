[Case]
  version = 3
  problem = steady
  geometry = axisymmetric_rz
[]

[Mesh]
  type = exodus
  file = ../moose/m33_multi_contact_rz_mesh.e
[]

[Materials]
  [solid]
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
        thermal_expansion = 0
        reference_temperature = 300
      []
    []
  []
[]

[Regions]
  [pellet]
    block = pellet
    material = solid
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [inner_clad]
    block = inner_clad
    material = solid
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [outer_clad]
    block = outer_clad
    material = solid
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]

[Contact]
  [pellet_to_inner_clad]
    primary = inner_clad_left
    secondary = pellet_right
    [mechanical]
      formulation = penalty
      penalty = 1e14
    []
  []
  [inner_to_outer_clad]
    primary = outer_clad_left
    secondary = inner_clad_right
    [mechanical]
      formulation = penalty
      penalty = 1e14
    []
  []
[]

[BoundaryConditions]
  [pellet_axis]
    type = dirichlet
    boundary = pellet_left
    field = radial_displacement
    value = 0
  []
  [inner_clad_inner]
    type = dirichlet
    boundary = inner_clad_left
    field = radial_displacement
    value = 0
  []
  [outer_clad_inner]
    type = dirichlet
    boundary = outer_clad_left
    field = radial_displacement
    value = 0
  []
  [pellet_closure]
    type = dirichlet
    boundary = pellet_right
    field = radial_displacement
    value = 1.1e-5
    scale_with_load = true
  []
  [inner_clad_closure]
    type = dirichlet
    boundary = inner_clad_right
    field = radial_displacement
    value = 1.1e-5
    scale_with_load = true
  []
  [pellet_bottom]
    type = dirichlet
    boundary = pellet_bottom
    field = axial_displacement
    value = 0
  []
  [inner_clad_bottom]
    type = dirichlet
    boundary = inner_clad_bottom
    field = axial_displacement
    value = 0
  []
  [outer_clad_bottom]
    type = dirichlet
    boundary = outer_clad_bottom
    field = axial_displacement
    value = 0
  []
  [pellet_top]
    type = dirichlet
    boundary = pellet_top
    field = axial_displacement
    value = 1e-5
  []
  [inner_clad_top]
    type = dirichlet
    boundary = inner_clad_top
    field = axial_displacement
    value = 1e-5
  []
  [outer_clad_top]
    type = dirichlet
    boundary = outer_clad_top
    field = axial_displacement
    value = 1e-5
  []
  [pellet_temperature]
    type = dirichlet
    boundary = pellet_left
    field = temperature
    value = 300
  []
  [inner_clad_temperature]
    type = dirichlet
    boundary = inner_clad_left
    field = temperature
    value = 300
  []
  [outer_clad_temperature]
    type = dirichlet
    boundary = outer_clad_left
    field = temperature
    value = 300
  []
[]

[Executioner]
  type = steady
  load_steps = 4
  cutback_factor = 0.5
  maximum_cutbacks = 12
  minimum_load_increment = 1e-6
[]
