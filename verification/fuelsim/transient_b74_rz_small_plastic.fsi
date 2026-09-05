[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/b7_rz_material.e
[]
[TimeFunctions]
  [extension]
    type = piecewise_linear
    times = 0 0.4 0.6 1
    values = 0 8e-6 7e-6 1e-5
  []
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
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 1e6
      hardening_modulus = 2e10
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
    type = dirichlet
    boundary = top
    field = axial_displacement
    value = 1
    function = extension
  []
[]
[Executioner]
  type = transient
  end_time = 1
  initial_time_step = 0.2
  minimum_time_step = 0.2
  maximum_time_step = 0.2
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
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
  csv = transient_b74_rz_small_plastic_summary.csv
  exodus = transient_b74_rz_small_plastic_results.e
[]
