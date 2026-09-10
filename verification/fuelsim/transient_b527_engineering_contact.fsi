[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = ../meshes/b527_engineering_contact.e
[]
[Materials]
  [fuel]
    [thermal]
      function = constant_thermophysical
      conductivity = 3
      density = 10000
      specific_heat = 300
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.30
    []
    [eigenstrains]
      [thermal]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 600
      []
    []
  []
  [clad]
    [thermal]
      function = constant_thermophysical
      conductivity = 15
      density = 6500
      specific_heat = 330
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e11
      poisson_ratio = 0.32
    []
    [eigenstrains]
      [thermal]
        function = isotropic_thermal_expansion
        thermal_expansion = 5e-6
        reference_temperature = 600
      []
    []
  []
[]
[Regions]
  [fuel]
    element = c3d8t
    block = fuel
    material = fuel
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 0
  []
  [clad]
    element = c3d8t
    block = clad
    material = clad
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[TimeFunctions]
  [fuel_temperature]
    type = piecewise_linear
    times = 0 10000
    values = 600 1200
  []
[]
[Contact]
  [fuel_clad_contact]
    primary = fuel_outer
    secondary = clad_inner
    [thermal]
      law = affine
      conductance = 100
      pressure_derivative = 1e-5
    []
    [mechanical]
      formulation = penalty
      penalty = 1e14
      mu = 0.1
      slip_tolerance = 0.005
      discretization = surface_to_surface
      sliding = finite
    []
  []
[]
[BoundaryConditions]
  [fuel_temperature]
    type = dirichlet
    boundary = fuel_inner
    field = temperature
    value = 1
    function = fuel_temperature
  []
  [coolant_temperature]
    type = dirichlet
    boundary = clad_outer
    field = temperature
    value = 600
  []
  [fuel_symmetry_y]
    type = dirichlet
    boundary = fuel_symmetry_y
    field = displacement_y
    value = 0
  []
  [fuel_symmetry_x]
    type = dirichlet
    boundary = fuel_symmetry_x
    field = displacement_x
    value = 0
  []
  [fuel_bottom]
    type = dirichlet
    boundary = fuel_bottom
    field = displacement_z
    value = 0
  []
  [clad_symmetry_y]
    type = dirichlet
    boundary = clad_symmetry_y
    field = displacement_y
    value = 0
  []
  [clad_symmetry_x]
    type = dirichlet
    boundary = clad_symmetry_x
    field = displacement_x
    value = 0
  []
  [clad_bottom]
    type = dirichlet
    boundary = clad_bottom
    field = displacement_z
    value = 0
  []
[]
[Executioner]
  type = transient
  end_time = 10000
  initial_time_step = 1000
  minimum_time_step = 1000
  maximum_time_step = 1000
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-11
  maximum_iterations = 60
  linear_solver = direct
  preconditioner = lu
  field_residual_scaling = true
  temperature_residual_absolute_tolerance = 1e-5
  mechanical_residual_absolute_tolerance = 1e-3
[]
[Outputs]
  console = false
  csv = transient_b527_engineering_contact_summary.csv
  exodus = transient_b527_engineering_contact_results.e
[]
