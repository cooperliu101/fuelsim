[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]

[Mesh]
  type = exodus
  file = ../moose/m31_transient_table_convection_rz_mesh.e
[]

[TimeFunctions]
  [power]
    type = piecewise_linear
    times = 0 2.5 5 10
    values = 0 1 0.5 1
  []
[]

[Materials]
  [fuel]
    [thermal]
      function = constant_thermophysical
      conductivity = 2
      density = 10000
      specific_heat = 300
    []

    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.3
    []
  []

[]

[Regions]
  [fuel]
    block_id = 0
    material = fuel
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 3e6
    heat_source_function = power
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
  [coolant]
    type = convection
    boundary = right
    heat_transfer_coefficient = 1000
    ambient_temperature = 500
  []
[]

[Executioner]
  type = transient
  end_time = 10
  initial_time_step = 3
  minimum_time_step = 0.125
  maximum_time_step = 3
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 3
  load_ramp_time = 0
[]

[Solver]
  absolute_tolerance = 1e-12
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 40
[]

[Outputs]
  console = false
[]
