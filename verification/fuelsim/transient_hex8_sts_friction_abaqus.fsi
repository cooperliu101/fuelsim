[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../abaqus/b40_hex8_sts_friction_mesh.e
[]

[Materials]
  [elastic]
    [thermal]
      function = constant_thermophysical
      conductivity = 1
      density = 1
      specific_heat = 1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0
    []
  []
[]

[Regions]
  [primary]
    element = c3d8t
    block = primary
    material = elastic
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary]
    element = c3d8t
    block = secondary
    material = elastic
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]

[TimeFunctions]
  [normal_path]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 -1e-4 -1e-4 -1e-4 -1e-4
  []
  [tangent_y_path]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 3e-6 12e-6 20e-6 25e-6
  []
  [tangent_z_path]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 4e-6 16e-6 15e-6 30e-6
  []
[]

[Contact]
  [interface]
    primary = primary_contact
    secondary = secondary_contact
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      sliding = small
      penalty = 1e8
      mu = 0.3
      slip_tolerance = 1e-5
    []
  []
[]

[BoundaryConditions]
  [primary_temperature]
    type = dirichlet
    boundary = primary_all
    field = temperature
    value = 300
  []
  [secondary_temperature]
    type = dirichlet
    boundary = secondary_all
    field = temperature
    value = 300
  []
  [primary_x]
    type = dirichlet
    boundary = primary_all
    field = displacement_x
    value = 0
  []
  [primary_y]
    type = dirichlet
    boundary = primary_all
    field = displacement_y
    value = 0
  []
  [primary_z]
    type = dirichlet
    boundary = primary_all
    field = displacement_z
    value = 0
  []
  [secondary_x]
    type = dirichlet
    boundary = secondary_all
    field = displacement_x
    value = 1
    function = normal_path
  []
  [secondary_y]
    type = dirichlet
    boundary = secondary_all
    field = displacement_y
    value = 1
    function = tangent_y_path
  []
  [secondary_z]
    type = dirichlet
    boundary = secondary_all
    field = displacement_z
    value = 1
    function = tangent_z_path
  []
[]

[Executioner]
  type = transient
  end_time = 4
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
  csv = transient_hex8_sts_friction_abaqus_summary.csv
  exodus = transient_hex8_sts_friction_abaqus_results.e
  checkpoint = transient_hex8_sts_friction_abaqus.checkpoint
  checkpoint_interval = 1
[]
