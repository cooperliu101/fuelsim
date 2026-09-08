[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = ../meshes/c3d20rt_probe.e
[]
[TimeFunctions]
  [pull]
    type = piecewise_linear
    times = 0 0.1 1
    values = 0 0.0015 0.006
  []
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 1
      density = 1
      specific_heat = 1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.3
    []
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 2e8
      hardening_modulus = 2e9
    []
  []
[]
[Regions]
  [solid]
    block = solid
    material = solid
    element = c3d20rt
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[BoundaryConditions]
  [temperature]
    type = dirichlet
    boundary = all
    field = temperature
    value = 600
  []
  [xzero]
    type = dirichlet
    boundary = xzero
    field = displacement_x
    value = 0
  []
  [yzero]
    type = dirichlet
    boundary = yzero
    field = displacement_y
    value = 0
  []
  [zzero]
    type = dirichlet
    boundary = zzero
    field = displacement_z
    value = 0
  []
  [pull]
    type = dirichlet
    boundary = xone
    field = displacement_x
    value = 1
    function = pull
  []
[]
[Executioner]
  type = transient
  end_time = 1
  initial_time_step = 0.1
  minimum_time_step = 0.1
  maximum_time_step = 0.1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-12
  maximum_iterations = 40
  field_residual_scaling = true
[]
[Outputs]
  console = true
  exodus = transient_c3d20rt_small_plastic_results.e
[]
