[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]

[Mesh]
  type = exodus
  file = ../moose/m21_transient_heat_rz_mesh.e
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
[]

[Executioner]
  type = transient
  end_time = 10
  initial_time_step = 1
  minimum_time_step = 1
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
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
