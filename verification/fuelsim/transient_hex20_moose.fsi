[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../moose/h20_06_transient.e
[]

[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 6000
      specific_heat = 1000
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.25
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 300
      []
    []
  []
[]

[Regions]
  [solid]
    element = c3d20t
    block_id = 0
    material = solid
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 6e6
  []
[]

[BoundaryConditions]
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
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-12
  step_tolerance = 1e-12
  maximum_iterations = 30
  linear_solver = direct
  preconditioner = lu
[]

[Outputs]
  console = false
  csv = transient_hex20_moose_summary.csv
  exodus = transient_hex20_moose_results.e
[]
