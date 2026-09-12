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
  [temperature]
    type = piecewise_linear
    times = 0 1
    values = 600 650
  []
  [force]
    type = piecewise_linear
    times = 0 1
    values = 0 28.274333882308138
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
  [temperature]
    type = dirichlet
    boundary = all_radial
    field = temperature
    value = 1
    function = temperature
  []
  [bottom]
    type = dirichlet
    boundary = bottom
    field = axial_displacement
    value = 0
  []
  [force]
    type = axial_force
    boundary = top
    value = 1
    function = force
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
  csv = transient_gps_uniform_finite_summary.csv
  exodus = transient_gps_uniform_finite_results.e
  history = transient_gps_uniform_finite_history.csv
  checkpoint = transient_gps_uniform_finite.checkpoint
[]
