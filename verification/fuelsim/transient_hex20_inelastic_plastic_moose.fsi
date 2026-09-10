[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = ../moose/h20_07_inelastic.e
[]
[Materials]
  [solid]
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
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 2e8
      hardening_modulus = 2e9
    []
  []
[]
[Regions]
  [solid]
    element = c3d20t
    block_id = 0
    material = solid
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[BoundaryConditions]
  [temperature]
    type = dirichlet
    boundary = left
    field = temperature
    value = 600
  []
  [fix_x]
    type = dirichlet
    boundary = left
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = bottom
    field = displacement_y
    value = 0
  []
  [fix_z]
    type = dirichlet
    boundary = back
    field = displacement_z
    value = 0
  []
  [pull_x]
    type = dirichlet
    boundary = right
    field = displacement_x
    value = 0.004
    scale_with_load = true
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
  load_ramp_time = 1
[]
[Solver]
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 40
  linear_solver = direct
  preconditioner = lu
[]
[Outputs]
  console = false
  csv = transient_hex20_inelastic_plastic_moose_summary.csv
  exodus = transient_hex20_inelastic_plastic_moose_results.e
[]
