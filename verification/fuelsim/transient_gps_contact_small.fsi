[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_1d
[]
[Mesh]
  type = exodus
  file = ../meshes/gps_two_slice.e
[]
[TimeFunctions]
  [closure]
    type = piecewise_linear
    times = 0 1 5 6 7 8 10
    values = 0 2e-5 2e-5 0 0 2e-5 2e-5
  []
  [translation]
    type = piecewise_linear
    times = 0 1 2 3 4 5 6 7 8 9 10
    values = 0 0 1e-8 2e-6 1.999e-6 -2e-6 -2e-6 2e-5 2e-5 2.001e-5 2.4e-5
  []
  [hot_temperature]
    type = piecewise_linear
    times = 0 1 10
    values = 600 700 700
  []
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 1000
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.3
    []
    [eigenstrains]
      [thermal]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 600
      []
    []
  []
[]
[Regions]
  [inner]
    block = inner
    element = cax2t_gps
    material = solid
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
  [outer]
    block = outer
    element = cax2t_gps
    material = solid
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[Contact]
  [layer1]
    primary = outer_interface1
    secondary = inner_interface1
    [thermal]
      law = affine
      conductance = 1e4
    []
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      sliding = small
      penalty = 1e12
      mu = 0.3
      slip_tolerance = 1e-5
    []
  []
  [layer2]
    primary = outer_interface2
    secondary = inner_interface2
    [thermal]
      law = affine
      conductance = 1e4
    []
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      sliding = small
      penalty = 1e12
      mu = 0.3
      slip_tolerance = 1e-5
    []
  []
[]
[BoundaryConditions]
  [hot]
    type = dirichlet
    boundary = inner_inside
    field = temperature
    value = 1
    function = hot_temperature
  []
  [cold]
    type = dirichlet
    boundary = outer_outside
    field = temperature
    value = 600
  []
  [inner_radial]
    type = dirichlet
    boundary = inner_radial
    field = radial_displacement
    value = 1
    function = closure
  []
  [outer_radial]
    type = dirichlet
    boundary = outer_radial
    field = radial_displacement
    value = 0
  []
  [inner_axial]
    type = dirichlet
    boundary = inner_axial
    field = axial_displacement
    value = 1
    function = translation
  []
  [outer_axial]
    type = dirichlet
    boundary = outer_axial
    field = axial_displacement
    value = 0
  []
[]
[Executioner]
  type = transient
  end_time = 10
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
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 40
[]
[Outputs]
  console = true
  csv = transient_gps_contact_small_summary.csv
  exodus = transient_gps_contact_small_results.e
  history = transient_gps_contact_small_history.csv
  checkpoint = transient_gps_contact_small.checkpoint
[]
