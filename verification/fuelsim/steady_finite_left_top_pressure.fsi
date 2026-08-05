[Case]
  version = 1
  problem = steady
[]

[Mesh]
  type = exodus
  file = ../moose/m42_left_top_pressure_rz_mesh.e
[]

[Regions]
  [solid]
    block_id = 0
    strain = finite
    conductivity_inverse_temperature = 0
    conductivity_constant = 1
    young_modulus = 2e11
    poisson_ratio = 0.3
    thermal_expansion = 0
    reference_temperature = 600
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]

[Contact]
[]

[BoundaryConditions]
  [right_radial]
    type = dirichlet
    boundary = right
    field = radial_displacement
    value = 0
  []
  [bottom_axial]
    type = dirichlet
    boundary = bottom
    field = axial_displacement
    value = 0
  []
  [temperature]
    type = dirichlet
    boundary = bottom
    field = temperature
    value = 600
  []
  [inner_pressure]
    type = pressure
    boundary = left
    value = 1e8
  []
  [top_pressure]
    type = pressure
    boundary = top
    value = 1e8
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
  maximum_iterations = 80
[]

[Outputs]
  console = false
[]
