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
  [inner_temperature]
    type = piecewise_linear
    times = 0 0.5 1
    values = 600 650 700
  []
  [outer_temperature]
    type = piecewise_linear
    times = 0 0.5 1
    values = 600 620 640
  []
  [inner_radial]
    type = piecewise_linear
    times = 0 0.5 1
    values = 0 0.0001 0.00015
  []
  [outer_radial]
    type = piecewise_linear
    times = 0 0.5 1
    values = 0 0.00005 0.00025
  []
  [top_axial]
    type = piecewise_linear
    times = 0 0.5 1
    values = 0 0.0002 -0.0001
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
  [solid]
    block = solid
    element = cax2t_gps
    material = solid
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[BoundaryConditions]
  [inner_temperature]
    type = dirichlet
    boundary = inner
    field = temperature
    value = 1
    function = inner_temperature
  []
  [outer_temperature]
    type = dirichlet
    boundary = outer
    field = temperature
    value = 1
    function = outer_temperature
  []
  [inner_radial]
    type = dirichlet
    boundary = inner
    field = radial_displacement
    value = 1
    function = inner_radial
  []
  [outer_radial]
    type = dirichlet
    boundary = outer
    field = radial_displacement
    value = 1
    function = outer_radial
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
    value = 1
    function = top_axial
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
  csv = transient_gps_nonuniform_finite_summary.csv
  exodus = transient_gps_nonuniform_finite_results.e
  history = transient_gps_nonuniform_finite_history.csv
  checkpoint = transient_gps_nonuniform_finite.checkpoint
[]
