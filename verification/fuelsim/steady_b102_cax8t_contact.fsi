[Case]
  version = 3
  problem = steady
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/b102_rz_contact.e
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
  [lower]
    block = lower
    material = solid
    element = cax8t
    strain = small
    initial_temperature = 650
    volumetric_heat_source = 0
  []
  [upper]
    block = upper
    material = solid
    element = cax8t
    strain = small
    initial_temperature = 650
    volumetric_heat_source = 0
  []
[]
[Contact]
  [interface]
    primary = upper_contact
    secondary = lower_contact
    [thermal]
      gap_conductivity = 0.1
      minimum_gap = 1e-6
    []
    [mechanical]
      formulation = penalty
      penalty = 1e13
    []
  []
[]
[BoundaryConditions]
  [lower_r]
    type = dirichlet
    boundary = lower_fixed
    field = radial_displacement
    value = 0
  []
  [upper_r]
    type = dirichlet
    boundary = upper_fixed
    field = radial_displacement
    value = 0
  []
  [lower_z]
    type = dirichlet
    boundary = lower_fixed
    field = axial_displacement
    value = 0
  []
  [upper_z]
    type = dirichlet
    boundary = upper_fixed
    field = axial_displacement
    value = -5e-6
  []
  [lower_t]
    type = dirichlet
    boundary = lower_temp
    field = temperature
    value = 700
  []
  [upper_t]
    type = dirichlet
    boundary = upper_temp
    field = temperature
    value = 600
  []
[]
[Executioner]
  type = steady
  load_steps = 1
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 50
[]
[Outputs]
  console = true
  csv = steady_b102_cax8t_contact_summary.csv
  exodus = steady_b102_cax8t_contact_results.e
[]
