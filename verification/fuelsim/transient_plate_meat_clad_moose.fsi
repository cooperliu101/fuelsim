[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../moose/b36_plate_meat_clad_mesh.e
[]

[Materials]
  [meat_material]
    [thermal]
      function = constant_thermophysical
      conductivity = 3
      density = 10970
      specific_heat = 300
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.316
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 600
      []
    []
    [creep]
      function = norton
      coefficient = 1e-7
      reference_stress = 1e8
      stress_exponent = 3
    []
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 2e8
      hardening_modulus = 2e9
    []
  []
  [clad_material]
    [thermal]
      function = constant_thermophysical
      conductivity = 16
      density = 6500
      specific_heat = 330
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.316
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 600
      []
    []
    [creep]
      function = norton
      coefficient = 1e-7
      reference_stress = 1e8
      stress_exponent = 3
    []
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 2e8
      hardening_modulus = 2e9
    []
  []
[]

[Regions]
  [meat]
    block = meat
    material = meat_material
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 3.291e5
  []
  [clad]
    block = clad
    material = clad_material
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 2.145e5
  []
[]

[TimeFunctions]
  [pull]
    type = piecewise_linear
    times = 0 0.5
    values = 0 2.4e-5
  []
[]

[BoundaryConditions]
  [fix_x]
    type = dirichlet
    boundary = plate_left
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = plate_bottom
    field = displacement_y
    value = 0
  []
  [fix_z]
    type = dirichlet
    boundary = plate_back
    field = displacement_z
    value = 0
  []
  [pull_x]
    type = dirichlet
    boundary = plate_right
    field = displacement_x
    value = 1
    function = pull
  []
[]

[Executioner]
  type = transient
  end_time = 0.5
  initial_time_step = 0.0125
  minimum_time_step = 0.0125
  maximum_time_step = 0.0125
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]

[Solver]
  linear_solver = direct
  preconditioner = lu
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 30
[]

[Outputs]
  console = false
  csv = transient_plate_meat_clad_moose_summary.csv
  exodus = transient_plate_meat_clad_moose_results.e
[]
