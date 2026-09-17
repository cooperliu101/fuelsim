[Case]
  version = 3
  problem = transient
  geometry = generalized_plane_strain
[]
[Mesh]
  type = exodus
  file = inelastic.e
[]
[TimeFunctions]
  [extension]
    type = piecewise_linear
    times = 0 1 2 3
    values = 0 0.001 0.001 0.0002
  []
[]
[Materials]
  [plastic]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e6
      poisson_ratio = 0.25
    []
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 1000
      hardening_modulus = 10000
    []
  []
  [creep]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e6
      poisson_ratio = 0.25
    []
    [creep]
      function = norton
      coefficient = 0.0001
      reference_stress = 1000
      stress_exponent = 1
    []
  []
  [coupled]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e6
      poisson_ratio = 0.25
    []
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 1000
      hardening_modulus = 10000
    []
    [creep]
      function = norton
      coefficient = 0.0001
      reference_stress = 1000
      stress_exponent = 1
    []
  []
[]
[Regions]
  [plastic]
    element = cpeg8t
    block = plastic
    material = plastic
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [creep]
    element = cpeg8t
    block = creep
    material = creep
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
  [coupled]
    element = cpeg8t
    block = coupled
    material = coupled
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]
[GeneralizedPlaneStrain]
  [plastic]
    blocks = plastic
    initial_thickness = 0.1
    u3 = 1
    u3_function = extension
    rotation_x = 0
    rotation_y = 0
  []
  [creep]
    blocks = creep
    initial_thickness = 0.1
    u3 = 1
    u3_function = extension
    rotation_x = 0
    rotation_y = 0
  []
  [coupled]
    blocks = coupled
    initial_thickness = 0.1
    u3 = 1
    u3_function = extension
    rotation_x = 0
    rotation_y = 0
  []
[]
[BoundaryConditions]
  [fix_x]
    type = dirichlet
    boundary = field
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = field
    field = displacement_y
    value = 0
  []
  [temperature]
    type = dirichlet
    boundary = corners
    field = temperature
    value = 300
  []
[]
[Executioner]
  type = transient
  end_time = 3
  initial_time_step = 0.25
  minimum_time_step = 0.25
  maximum_time_step = 0.25
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
  include_thermal_time_term = false
[]
[Solver]
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-11
  step_tolerance = 1e-14
  maximum_iterations = 30
  linear_solver = direct
  direct_factorization = mumps
[]
[Outputs]
  console = true
  exodus = inelastic_history_results.e
  history = inelastic_history_history.csv
  checkpoint = inelastic_history.chk
[]
