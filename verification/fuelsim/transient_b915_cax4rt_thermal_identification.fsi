[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/b915_rz_thermal.e
[]
[TimeFunctions]
  [temperature2]
    type = piecewise_linear
    times = 0 0.1 0.2
    values = 600 620 620
  []
  [temperature3]
    type = piecewise_linear
    times = 0 0.1 0.2
    values = 600 590 610
  []
  [temperature4]
    type = piecewise_linear
    times = 0 0.1 0.2
    values = 600 610 590
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
      young_modulus = 1e6
      poisson_ratio = 0.25
    []
  []
[]
[Regions]
  [solid]
    block = solid
    material = solid
    element = cax4rt
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[BoundaryConditions]
  [radial]
    type = dirichlet
    boundary = all
    field = radial_displacement
    value = 0
  []
  [axial]
    type = dirichlet
    boundary = all
    field = axial_displacement
    value = 0
  []
  [t1]
    type = dirichlet
    boundary = t1
    field = temperature
    value = 600
  []
  [t2]
    type = dirichlet
    boundary = t2
    field = temperature
    value = 1
    function = temperature2
  []
  [t3]
    type = dirichlet
    boundary = t3
    field = temperature
    value = 1
    function = temperature3
  []
  [t4]
    type = dirichlet
    boundary = t4
    field = temperature
    value = 1
    function = temperature4
  []
[]
[Executioner]
  type = transient
  include_thermal_time_term = false
  end_time = 0.2
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
  maximum_iterations = 40
[]
[Outputs]
  console = true
  csv = transient_b915_cax4rt_thermal_identification_summary.csv
  exodus = transient_b915_cax4rt_thermal_identification_results.e
[]
