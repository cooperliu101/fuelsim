[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/cax4t_material_temperature.e
[]
[TimeFunctions]
  [extension]
    type = piecewise_linear
    times = 0 0.1 1
    values = 0 0.0001 0.0001
  []
  [temperature1]
    type = piecewise_linear
    times = 0 0.1 1
    values = 600 600 600
  []
  [temperature2]
    type = piecewise_linear
    times = 0 0.1 1
    values = 600 650 650
  []
  [temperature3]
    type = piecewise_linear
    times = 0 0.1 1
    values = 600 700 700
  []
  [temperature4]
    type = piecewise_linear
    times = 0 0.1 1
    values = 600 750 750
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
    [creep]
      function = linear_temperature_norton
      coefficient = 1e-19
      reference_stress = 1
      stress_exponent = 2.2
      reference_temperature = 600
      coefficient_temperature_coefficient = 0
      reference_stress_temperature_coefficient = 0
      stress_exponent_temperature_coefficient = 0.001
    []
  []
[]
[Regions]
  [solid]
    block = solid
    element = cax4t
    material = solid
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[BoundaryConditions]
  [r1]
    type = dirichlet
    boundary = n1
    field = radial_displacement
    value = 0
  []
  [r2]
    type = dirichlet
    boundary = n2
    field = radial_displacement
    value = 0
  []
  [r3]
    type = dirichlet
    boundary = n3
    field = radial_displacement
    value = 0
  []
  [r4]
    type = dirichlet
    boundary = n4
    field = radial_displacement
    value = 0
  []
  [z1]
    type = dirichlet
    boundary = n1
    field = axial_displacement
    value = 0
  []
  [z2]
    type = dirichlet
    boundary = n2
    field = axial_displacement
    value = 0
  []
  [z3]
    type = dirichlet
    boundary = n3
    field = axial_displacement
    value = 1
    function = extension
  []
  [z4]
    type = dirichlet
    boundary = n4
    field = axial_displacement
    value = 1
    function = extension
  []
  [t1]
    type = dirichlet
    boundary = n1
    field = temperature
    value = 1
    function = temperature1
  []
  [t2]
    type = dirichlet
    boundary = n2
    field = temperature
    value = 1
    function = temperature2
  []
  [t3]
    type = dirichlet
    boundary = n3
    field = temperature
    value = 1
    function = temperature3
  []
  [t4]
    type = dirichlet
    boundary = n4
    field = temperature
    value = 1
    function = temperature4
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
  csv = transient_cax4t_creep_exponent_finite_summary.csv
  exodus = transient_cax4t_creep_exponent_finite_results.e
[]
