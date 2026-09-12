[Case]
  version = 3
  problem = steady
  geometry = axisymmetric_1d
[]
[Mesh]
  type = exodus
  file = ../meshes/gps_uniform.e
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 1000
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e9
      poisson_ratio = 0.3
    []
    [eigenstrains]
      [thermal]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-4
        reference_temperature = 600
      []
    []
  []
[]
[Regions]
  [solid]
    block = solid
    element = cax2t_gps
    material = solid
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[BoundaryConditions]
  [temperature]
    type = dirichlet
    boundary = all_radial
    field = temperature
    value = 1100
  []
  [bottom]
    type = dirichlet
    boundary = bottom
    field = axial_displacement
    value = 0
  []
[]
[Executioner]
  type = steady
  load_steps = 5
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
  csv = steady_gps_uniform_finite_summary.csv
  exodus = steady_gps_uniform_finite_results.e
[]
