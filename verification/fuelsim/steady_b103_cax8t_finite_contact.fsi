[Case]
  version = 3
  problem = steady
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/b10_cax8t_contact.e
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
    element = cax8t
    material = solid
    strain = finite
    initial_temperature = 650
    volumetric_heat_source = 0
  []
  [upper]
    block = upper
    element = cax8t
    material = solid
    strain = finite
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
  [radial]
    type = dirichlet
    boundary = lower_radial
    field = radial_displacement
    value = 0
  []
  [upper_radial]
    type = dirichlet
    boundary = upper_radial
    field = radial_displacement
    value = 0
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
    value = -5e-6
  []
  [hot]
    type = dirichlet
    boundary = bottom
    field = temperature
    value = 700
  []
  [cold]
    type = dirichlet
    boundary = top
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
  csv = steady_b103_cax8t_finite_contact_summary.csv
  exodus = steady_b103_cax8t_finite_contact_results.e
[]
