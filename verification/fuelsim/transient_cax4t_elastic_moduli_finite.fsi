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
  [inner_radial]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.00004 0.00012
  []
  [outer_radial]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.00005 0.00015
  []
  [top_axial]
    type = piecewise_linear
    times = 0 1 2
    values = 0 0.00025 0.00045
  []
  [temperature1]
    type = piecewise_linear
    times = 0 1 1.25 2
    values = 600 600 750 750
  []
  [temperature2]
    type = piecewise_linear
    times = 0 1 1.25 2
    values = 650 650 700 700
  []
  [temperature3]
    type = piecewise_linear
    times = 0 1 1.25 2
    values = 700 700 650 650
  []
  [temperature4]
    type = piecewise_linear
    times = 0 1 1.25 2
    values = 750 750 600 600
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
      function = linear_temperature_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.2
      reference_temperature = 600
      young_modulus_temperature_coefficient = -2e6
      poisson_ratio_temperature_coefficient = 0.001
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
    value = 1
    function = inner_radial
  []
  [r2]
    type = dirichlet
    boundary = n2
    field = radial_displacement
    value = 1
    function = outer_radial
  []
  [r3]
    type = dirichlet
    boundary = n3
    field = radial_displacement
    value = 1
    function = outer_radial
  []
  [r4]
    type = dirichlet
    boundary = n4
    field = radial_displacement
    value = 1
    function = inner_radial
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
    function = top_axial
  []
  [z4]
    type = dirichlet
    boundary = n4
    field = axial_displacement
    value = 1
    function = top_axial
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
  initial_time_step = 0.25
  minimum_time_step = 0.25
  maximum_time_step = 0.25
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
  csv = transient_cax4t_elastic_moduli_finite_summary.csv
  exodus = transient_cax4t_elastic_moduli_finite_results.e
[]
