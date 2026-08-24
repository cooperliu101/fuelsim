[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../abaqus/b43_hex8_finite_strain_contact_mesh.e
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
      young_modulus = 1e7
      poisson_ratio = 0
    []
  []
[]

[Regions]
  [primary]
    block = primary
    material = elastic
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [secondary]
    block = secondary
    material = elastic
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]

[TimeFunctions]
  [primary_y0_x]
    type = piecewise_linear
    times = 0 20
    values = 0 -0.015
  []
  [primary_y1_x]
    type = piecewise_linear
    times = 0 20
    values = 0 0.015
  []
  [primary_y0_y]
    type = piecewise_linear
    times = 0 20
    values = 0 -0.05
  []
  [primary_y1_y]
    type = piecewise_linear
    times = 0 20
    values = 0 0.05
  []
  [primary_z0_z]
    type = piecewise_linear
    times = 0 20
    values = 0 -0.06
  []
  [primary_z1_z]
    type = piecewise_linear
    times = 0 20
    values = 0 0.06
  []
  [secondary_y0_x]
    type = piecewise_linear
    times = 0 20
    values = 0 -0.095
  []
  [secondary_y1_x]
    type = piecewise_linear
    times = 0 20
    values = 0 -0.065
  []
  [secondary_y0_y]
    type = piecewise_linear
    times = 0 20
    values = 0 -0.05
  []
  [secondary_y1_y]
    type = piecewise_linear
    times = 0 20
    values = 0 0.05
  []
  [secondary_z0_z]
    type = piecewise_linear
    times = 0 20
    values = 0 -0.02
  []
  [secondary_z1_z]
    type = piecewise_linear
    times = 0 20
    values = 0 0.06
  []
[]

[Contact]
  [interface]
    primary = primary_contact
    secondary = secondary_contact
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      sliding = finite
      penalty = 1e8
      mu = 0
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
  [primary_y0_x]
    type = dirichlet
    boundary = primary_y0
    field = displacement_x
    value = 1
    function = primary_y0_x
  []
  [primary_y1_x]
    type = dirichlet
    boundary = primary_y1
    field = displacement_x
    value = 1
    function = primary_y1_x
  []
  [primary_y0_y]
    type = dirichlet
    boundary = primary_y0
    field = displacement_y
    value = 1
    function = primary_y0_y
  []
  [primary_y1_y]
    type = dirichlet
    boundary = primary_y1
    field = displacement_y
    value = 1
    function = primary_y1_y
  []
  [primary_z0_z]
    type = dirichlet
    boundary = primary_z0
    field = displacement_z
    value = 1
    function = primary_z0_z
  []
  [primary_z1_z]
    type = dirichlet
    boundary = primary_z1
    field = displacement_z
    value = 1
    function = primary_z1_z
  []
  [secondary_y0_x]
    type = dirichlet
    boundary = secondary_y0
    field = displacement_x
    value = 1
    function = secondary_y0_x
  []
  [secondary_y1_x]
    type = dirichlet
    boundary = secondary_y1
    field = displacement_x
    value = 1
    function = secondary_y1_x
  []
  [secondary_y0_y]
    type = dirichlet
    boundary = secondary_y0
    field = displacement_y
    value = 1
    function = secondary_y0_y
  []
  [secondary_y1_y]
    type = dirichlet
    boundary = secondary_y1
    field = displacement_y
    value = 1
    function = secondary_y1_y
  []
  [secondary_z0_z]
    type = dirichlet
    boundary = secondary_z0
    field = displacement_z
    value = 1
    function = secondary_z0_z
  []
  [secondary_z1_z]
    type = dirichlet
    boundary = secondary_z1
    field = displacement_z
    value = 1
    function = secondary_z1_z
  []
[]

[Executioner]
  type = transient
  end_time = 20
  initial_time_step = 1
  minimum_time_step = 1
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
  include_thermal_time_term = false
[]

[Solver]
  absolute_tolerance = 1e-7
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
