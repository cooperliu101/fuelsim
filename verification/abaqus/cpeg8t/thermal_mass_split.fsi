[Case]
  version = 3
  problem = transient
  geometry = generalized_plane_strain
[]
[Mesh]
  type = exodus
  file = distorted.e
[]
[TimeFunctions]
  [ramp]
    type = piecewise_linear
    times = 0 1
    values = 0 1
  []
[]
[Materials]
  [solid]
    [thermal]
      function = linear_temperature_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
      reference_temperature = 300
      conductivity_temperature_coefficient = 0.01
      specific_heat_temperature_coefficient = 0.1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e6
      poisson_ratio = 0.25
    []
  []
[]
[Regions]
  [solid]
    element = cpeg8t
    block = solid
    material = solid
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 1e6
    body_acceleration = 2 -9.81
  []
[]
[GeneralizedPlaneStrain]
  [solid]
    blocks = solid
    initial_thickness = 0.1
    u3 = 0.01
    u3_function = ramp
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
  [top_pressure]
    type = pressure
    boundary = top
    value = 100
    configuration = current
  []
  [right_traction]
    type = traction
    boundary = right
    field = displacement_x
    value = 50
    configuration = current
  []
  [bottom_traction]
    type = traction
    boundary = bottom
    field = displacement_y
    value = 30
    configuration = reference
  []
[]
[Executioner]
  type = transient
  end_time = 0.5
  initial_time_step = 0.25
  minimum_time_step = 0.25
  maximum_time_step = 0.25
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
  include_thermal_time_term = true
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
  exodus = thermal_mass_split_results.e
  history = thermal_mass_split_history.csv
  checkpoint = thermal_mass_split.chk
  checkpoint_interval = 1
[]
