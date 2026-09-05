[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/b7_rz_material.e
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.3
    []
    [creep]
      function = norton
      coefficient = 1e-4
      reference_stress = 1e8
      stress_exponent = 3
    []
  []
[]
[Regions]
  [solid]
    block_id = 0
    material = solid
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
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
  [top]
    type = traction
    boundary = top
    field = axial_displacement
    value = 1e8
    configuration = reference
    scale_with_load = true
  []
[]
[Executioner]
  type = transient
  end_time = 1.0000000001
  initial_time_step = 1e-10
  minimum_time_step = 1e-10
  maximum_time_step = 0.2
  growth_factor = 2e9
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 1e-10
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 40
[]
[Outputs]
  console = true
  csv = transient_b72_rz_small_creep_summary.csv
  exodus = transient_b72_rz_small_creep_results.e
[]
