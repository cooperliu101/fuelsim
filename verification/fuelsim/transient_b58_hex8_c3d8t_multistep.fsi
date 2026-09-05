[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = ../meshes/b55_hex8_two_elements.e
[]
[Materials]
  [solid]
    [thermal]
      function = linear_temperature_thermophysical
      conductivity = 4
      density = 2000
      specific_heat = 3000
      reference_temperature = 300
      conductivity_temperature_coefficient = 0.01
      density_temperature_coefficient = -1
      specific_heat_temperature_coefficient = 4
    []
    [elasticity]
      function = linear_temperature_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.25
      reference_temperature = 300
      young_modulus_temperature_coefficient = -1e8
      poisson_ratio_temperature_coefficient = 1e-4
    []
    [eigenstrains]
      [thermal_expansion]
        function = linear_temperature_isotropic_thermal_expansion
        thermal_expansion = 1.2e-5
        reference_temperature = 300
        thermal_expansion_temperature_coefficient = 2e-8
      []
    []
  []
[]
[Regions]
  [solid]
    block = solid
    material = solid
    strain = small
    element = c3d8t
    initial_temperature = 300
    volumetric_heat_source = 1
    heat_source_function = body_heat
  []
[]
[TimeFunctions]
  [right_temperature]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 300 360 420 420 330
  []
  [body_heat]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 2e6 4e6 0 0
  []
  [surface_heat]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 2000 5000 5000 0
  []
[]
[BoundaryConditions]
  [left_temperature]
    type = dirichlet
    boundary = x0
    field = temperature
    value = 300
  []
  [right_temperature]
    type = dirichlet
    boundary = x2
    field = temperature
    value = 1
    function = right_temperature
  []
  [fix_x]
    type = dirichlet
    boundary = x0
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = y0
    field = displacement_y
    value = 0
  []
  [fix_z]
    type = dirichlet
    boundary = z0
    field = displacement_z
    value = 0
  []
  [top_heat_flux]
    type = heat_flux
    boundary = z1
    value = 1
    function = surface_heat
  []
  [side_convection]
    type = convection
    boundary = y1
    heat_transfer_coefficient = 20
    ambient_temperature = 280
  []
[]
[Executioner]
  type = transient
  end_time = 4
  initial_time_step = 1
  minimum_time_step = 1
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-13
  step_tolerance = 1e-12
  maximum_iterations = 20
  linear_solver = direct
  preconditioner = lu
  line_search = basic
  backtracking_fallback = true
  field_residual_scaling = true
  residual_reduction_tolerance = 1e-11
  temperature_residual_absolute_tolerance = 1e-8
  mechanical_residual_absolute_tolerance = 1e-3
[]
[Outputs]
  console = true
  csv = transient_b58_hex8_c3d8t_multistep_summary.csv
  exodus = transient_b58_hex8_c3d8t_multistep_results.e
[]
