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
  [loading]
    type = piecewise_linear
    times = 0 1
    values = 0 1
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
  [interface]
    primary = outer_interface
    secondary = inner_interface
    [thermal]
      law = affine
      conductance = 1e4
    []
    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      sliding = small
      penalty = 1e10
      mu = 0.3
      slip_tolerance = 1e-6
    []
  []
[]
[BoundaryConditions]
  [inner_temperature]
    type = dirichlet
    boundary = inner_radial
    field = temperature
    value = 600
  []
  [outer_temperature]
    type = dirichlet
    boundary = outer_radial
    field = temperature
    value = 600
  []
  [inner_inside_radial1]
    type = dirichlet
    boundary = inner_inside1
    field = radial_displacement
    value = 16e-6
    function = loading
  []
  [inner_interface_radial1]
    type = dirichlet
    boundary = inner_interface1
    field = radial_displacement
    value = 20e-6
    function = loading
  []
  [inner_inside_radial2]
    type = dirichlet
    boundary = inner_inside2
    field = radial_displacement
    value = 15.2e-6
    function = loading
  []
  [inner_interface_radial2]
    type = dirichlet
    boundary = inner_interface2
    field = radial_displacement
    value = 19e-6
    function = loading
  []
  [outer_radial]
    type = dirichlet
    boundary = outer_radial
    field = radial_displacement
    value = 0
  []
  [inner_bottom]
    type = dirichlet
    boundary = inner_bottom
    field = axial_displacement
    value = 0
  []
  [inner_top]
    type = dirichlet
    boundary = inner_top
    field = axial_displacement
    value = 1e-4
    function = loading
  []
  [outer_bottom]
    type = dirichlet
    boundary = outer_bottom
    field = axial_displacement
    value = 0
  []
  [outer_top]
    type = dirichlet
    boundary = outer_top
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
  csv = transient_gps_chain_small_summary.csv
  exodus = transient_gps_chain_small_results.e
  history = transient_gps_chain_small_history.csv
  checkpoint = transient_gps_chain_small.checkpoint
[]
