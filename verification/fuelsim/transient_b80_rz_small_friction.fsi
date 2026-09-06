[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/b8_rz_sliding.e
[]
[TimeFunctions]
  [closure]
    type = piecewise_linear
    times = 0 0.1 1
    values = 0 5e-6 5e-6
  []
  [translation]
    type = piecewise_linear
    times = 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
    values = 0 0 1e-9 1e-5 2e-5 0.000019999 1e-5 0 -1e-5 -0.000009999 0
  []
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.3
    []
  []
[]
[Regions]
  [inner]
    block = inner
    element = cax4t
    material = solid
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
  [outer]
    block = outer
    element = cax4t
    material = solid
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[Contact]
  [interface]
    primary = outer_left
    secondary = inner_right
    [mechanical]
      formulation = penalty
      discretization = node_to_surface
      sliding = finite
      penalty = 1e13
      mu = 0.2
      slip_tolerance = 0.001
    []
  []
[]
[BoundaryConditions]
  [inner_radial]
    type = dirichlet
    boundary = inner_left
    field = radial_displacement
    value = 1
    function = closure
  []
  [inner_axial]
    type = dirichlet
    boundary = inner_left
    field = axial_displacement
    value = 1
    function = translation
  []
  [outer_radial]
    type = dirichlet
    boundary = outer_right
    field = radial_displacement
    value = 0
  []
  [outer_axial]
    type = dirichlet
    boundary = outer_right
    field = axial_displacement
    value = 0
  []
[]
[Executioner]
  type = transient
  end_time = 1
  initial_time_step = 0.1
  minimum_time_step = 0.1
  maximum_time_step = 0.1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-10
  step_tolerance = 1e-16
  maximum_iterations = 60
[]
[Outputs]
  console = true
  csv = transient_b80_rz_small_friction_summary.csv
  exodus = transient_b80_rz_small_friction_results.e
[]
