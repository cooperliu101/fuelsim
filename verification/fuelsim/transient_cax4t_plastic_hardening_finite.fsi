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
    times = 0 1 2
    values = 0 0.0004 0.0008
  []
  [temperature1]
    type = piecewise_linear
    times = 0 1 2
    values = 600 600 750
  []
  [temperature2]
    type = piecewise_linear
    times = 0 1 2
    values = 650 650 700
  []
  [temperature3]
    type = piecewise_linear
    times = 0 1 2
    values = 700 700 650
  []
  [temperature4]
    type = piecewise_linear
    times = 0 1 2
    values = 750 750 600
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
    [plasticity]
      function = linear_temperature_isotropic_hardening
      yield_stress = 2e6
      hardening_modulus = 2e7
      reference_temperature = 600
      yield_stress_temperature_coefficient = 0
      hardening_temperature_coefficient = 4e5
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
  end_time = 2
  initial_time_step = 0.2
  minimum_time_step = 0.2
  maximum_time_step = 0.2
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
  csv = transient_cax4t_plastic_hardening_finite_summary.csv
  exodus = transient_cax4t_plastic_hardening_finite_results.e
[]
