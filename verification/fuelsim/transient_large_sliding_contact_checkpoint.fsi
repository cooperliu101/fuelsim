[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../moose/m52_large_sliding_contact_rz_mesh.e
[]
[Materials]
  [lower]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 10000
      specific_heat = 300
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.3
    []
  []
  [upper]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 10000
      specific_heat = 300
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.3
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
  [lower_outer]
    type = dirichlet
    boundary = lower_right
    field = radial_displacement
    value = 0.00023
    scale_with_load = true
  []
  [upper_axis]
    type = dirichlet
    boundary = upper_left
    field = radial_displacement
    value = 0
  []
  [upper_outer]
    type = dirichlet
    boundary = upper_right
    field = radial_displacement
    value = 0.00004
    scale_with_load = true
  []
  [lower_bottom]
    type = dirichlet
    boundary = lower_bottom
    field = axial_displacement
    value = 0
  []
  [lower_contact_displacement]
    type = dirichlet
    boundary = lower_top
    field = axial_displacement
    value = 3e-6
    scale_with_load = true
  []
  [upper_bottom]
    type = dirichlet
    boundary = upper_bottom
    field = axial_displacement
    value = 0
  []
  [upper_top]
    type = dirichlet
    boundary = upper_top
    field = axial_displacement
    value = -1e-6
    scale_with_load = true
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
  type = transient
  end_time = 0.5
  initial_time_step = 0.05
  minimum_time_step = 0.05
  maximum_time_step = 0.05
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 2
  load_ramp_time = 1
[]
[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 80
[]
[Outputs]
  console = false
  checkpoint = transient_large_sliding_contact_split.checkpoint
[]
