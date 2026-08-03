[Case]
  version = 1
  problem = steady
[]

[Mesh]
  type = exodus
  file = ../moose/m0_simple_fuel_rz_mesh.e
[]

[Regions]
  [fuel]
    block_id = 0
    strain = small
    conductivity_inverse_temperature = 3824
    conductivity_constant = 0.61
    young_modulus = 2e11
    poisson_ratio = 0.316
    thermal_expansion = 1e-5
    reference_temperature = 600
    initial_temperature = 600
    volumetric_heat_source = 2e8
  []
[]

[Contact]
[]

[BoundaryConditions]
  [axis]
    type = dirichlet
    boundary = fuel_left
    field = radial_displacement
    value = 0
  []
  [bottom]
    type = dirichlet
    boundary = fuel_bottom
    field = axial_displacement
    value = 0
  []
  [outer_temperature]
    type = dirichlet
    boundary = fuel_right
    field = temperature
    value = 600
  []
[]

[Executioner]
  type = steady
  load_steps = 1
[]

[Solver]
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 40
[]

[Outputs]
  console = false
[]
