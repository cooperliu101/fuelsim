[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../abaqus/h20_30_hex20_sts_quadratic_cylinder_mesh.e
[]

[Materials]
  [primary]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1
      specific_heat = 1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.25
    []
  []
  [secondary]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1
      specific_heat = 1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.25
    []
  []
[]

[Regions]
  [primary]
    block = primary
    material = primary
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary]
    block = secondary
    material = secondary
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]

[TimeFunctions]
  [axial_displacement_path]
    type = piecewise_linear
    times = 0 1 2 3 4 5 6 7
    values = 0 2e-6 4e-6 12e-6 24e-6 4e-6 -20e-6 -18e-6
  []
[]

[Contact]
  [interface]
    primary = primary_contact
    secondary = secondary_contact
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      penalty = 1e11
      mu = 0.3
      slip_tolerance = 1e-5
    []
  []
[]

[BoundaryConditions]
  [primary_temperature]
    type = dirichlet
    boundary = primary_inner
    field = temperature
    value = 300
  []
  [secondary_temperature]
    type = dirichlet
    boundary = secondary_outer
    field = temperature
    value = 300
  []
  [primary_fix_x]
    type = dirichlet
    boundary = primary_inner
    field = displacement_x
    value = 0
  []
  [primary_fix_y]
    type = dirichlet
    boundary = primary_inner
    field = displacement_y
    value = 0
  []
  [primary_fix_z]
    type = dirichlet
    boundary = primary_back
    field = displacement_z
    value = 0
  []
  [secondary_symmetry_lower]
    type = dirichlet
    boundary = secondary_theta_lower
    field = displacement_y
    value = 0
  []
  [secondary_symmetry_upper]
    type = dirichlet
    boundary = secondary_theta_upper
    field = displacement_x
    value = 0
  []
  [outer_pressure]
    type = pressure
    boundary = secondary_outer
    value = 10000
    configuration = reference
  []
  [outer_axial_displacement]
    type = dirichlet
    boundary = secondary_outer
    field = displacement_z
    value = 1
    function = axial_displacement_path
  []
[]

[Executioner]
  type = transient
  end_time = 7
  initial_time_step = 1
  minimum_time_step = 0.125
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 3
  load_ramp_time = 0
[]

[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-11
  step_tolerance = 1e-12
  maximum_iterations = 40
  linear_solver = direct
  direct_factorization = mumps
  field_residual_scaling = true
  linear_relative_tolerance = 1e-11
  maximum_linear_iterations = 400
[]

[Outputs]
  console = true
[]
