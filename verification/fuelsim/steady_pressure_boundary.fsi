[Case]
  version = 3
  problem = steady
  geometry = axisymmetric_rz
[]

[Mesh]
  type = exodus
  file = ../moose/m22_coupled_plastic_creep_traction_rz_mesh.e
[]

[Materials]
  [material]
    [thermal]
      function = constant_thermophysical
      conductivity = 1
      density = 1
      specific_heat = 1
    []

    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.3
    []
  []

[]

[Regions]
  [material]
    block_id = 0
    material = material
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]

[BoundaryConditions]
  [axis]
    type = dirichlet
    boundary = left
    field = radial_displacement
    value = 0
  []
  [bottom]
    type = dirichlet
    boundary = bottom
    field = axial_displacement
    value = 0
  []
  [temperature]
    type = dirichlet
    boundary = top
    field = temperature
    value = 600
  []
  [outer_pressure]
    type = pressure
    boundary = right
    value = 1e6
    scale_with_load = true
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
  csv = steady_pressure_boundary_summary.csv
  exodus = steady_pressure_boundary_results.e
[]
