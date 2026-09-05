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
      function = constant_thermophysical
      conductivity = 4
      density = 2000
      specific_heat = 3000
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.25
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1.2e-5
        reference_temperature = 300
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
    volumetric_heat_source = 4e6
  []
[]
[TimeFunctions]
  [right_temperature]
    type = piecewise_linear
    times = 0 1
    values = 300 400
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
    value = 5000
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
  end_time = 1
  initial_time_step = 1
  minimum_time_step = 1
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 1
[]
[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-12
  step_tolerance = 1e-12
  maximum_iterations = 12
  linear_solver = direct
  preconditioner = lu
  line_search = basic
  backtracking_fallback = true
  field_residual_scaling = true
  residual_reduction_tolerance = 1e-9
  temperature_residual_absolute_tolerance = 1e-7
  mechanical_residual_absolute_tolerance = 1e-3
[]
[Outputs]
  console = true
  csv = transient_b55_hex8_c3d8t_thermal_summary.csv
  exodus = transient_b55_hex8_c3d8t_thermal_results.e
[]
