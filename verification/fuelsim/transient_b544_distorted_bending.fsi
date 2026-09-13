[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = ../meshes/b544_distorted_bending.e
[]
[TimeFunctions]
  [right_temperature]
    type = piecewise_linear
    times = 0 1000000
    values = 300 700
  []
[]
[Materials]
  [solid]
    [thermal]
      function = linear_temperature_thermophysical
      conductivity = 10
      density = 2000
      specific_heat = 500
      reference_temperature = 300
      conductivity_temperature_coefficient = 0.02
      specific_heat_temperature_coefficient = 1.25
    []
    [elasticity]
      function = linear_temperature_isotropic
      young_modulus = 2e8
      poisson_ratio = 0.25
      reference_temperature = 300
      young_modulus_temperature_coefficient = -2e5
      poisson_ratio_temperature_coefficient = 0
    []
    [eigenstrains]
      [thermal]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 300
      []
    []
  []
[]
[Regions]
  [solid]
    block = solid
    material = solid
    strain = finite
    element = c3d8rt
    initial_temperature = 300
    volumetric_heat_source = 0.5
  []
[]
[BoundaryConditions]
  [left_temperature]
    type = dirichlet
    boundary = left
    field = temperature
    value = 300
  []
  [right_temperature]
    type = dirichlet
    boundary = right
    field = temperature
    value = 1
    function = right_temperature
  []
  [fix_x]
    type = dirichlet
    boundary = left
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = left
    field = displacement_y
    value = 0
  []
  [fix_z]
    type = dirichlet
    boundary = left
    field = displacement_z
    value = 0
  []
  [top_heat_flux]
    type = heat_flux
    boundary = z_high
    value = 50
  []
  [tip_bending_traction]
    type = traction
    boundary = right
    field = displacement_z
    value = -2e4
  []
  [side_convection]
    type = convection
    boundary = y_high
    heat_transfer_coefficient = 5
    ambient_temperature = 280
  []
[]
[Executioner]
  type = transient
  end_time = 1000000
  initial_time_step = 100000
  minimum_time_step = 100000
  maximum_time_step = 100000
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-12
  maximum_iterations = 50
  linear_solver = direct
  preconditioner = lu
  field_residual_scaling = true
  residual_reduction_tolerance = 1e-10
  temperature_residual_absolute_tolerance = 1e-6
  mechanical_residual_absolute_tolerance = 1e-4
[]
[Outputs]
  console = false
  csv = transient_b544_distorted_bending_summary.csv
  exodus = transient_b544_distorted_bending_results.e
[]
