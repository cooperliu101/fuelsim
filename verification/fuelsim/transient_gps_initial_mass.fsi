[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_1d
[]
[Mesh]
  type = exodus
  file = ../meshes/gps_uniform.e
[]
[TimeFunctions]
  [ramp]
    type = piecewise_linear
    times = 0 1
    values = 0 1
  []
[]
[Materials]
  [solid]
    [thermal]
      function = linear_temperature_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 1000
      reference_temperature = 600
      conductivity_temperature_coefficient = 0.1
      density_temperature_coefficient = 2
      specific_heat_temperature_coefficient = 5
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.3
    []
  []
[]
[Regions]
  [solid]
    block = solid
    element = cax2t_gps
    material = solid
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 1e7
  []
[]
[BoundaryConditions]
  [inner]
    type = dirichlet
    boundary = inner
    field = radial_displacement
    value = 0.0008
    function = ramp
  []
  [outer]
    type = dirichlet
    boundary = outer
    field = radial_displacement
    value = 0.001
    function = ramp
  []
  [bottom]
    type = dirichlet
    boundary = bottom
    field = axial_displacement
    value = 0
  []
  [top]
    type = dirichlet
    boundary = top
    field = axial_displacement
    value = 0.001
    function = ramp
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
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 40
[]
[Outputs]
  console = true
  csv = transient_gps_initial_mass_summary.csv
  exodus = transient_gps_initial_mass_results.e
  history = transient_gps_initial_mass_history.csv
  checkpoint = transient_gps_initial_mass.checkpoint
[]
