[Case]
  version = 3
  problem = steady
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/b10_cax8t_material.e
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 2000
      specific_heat = 1000
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e6
      poisson_ratio = 0.25
    []
    [eigenstrains]
      [expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-3
        reference_temperature = 300
      []
    []
  []
[]
[Regions]
  [solid]
    block = solid
    material = solid
    element = cax8rt
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
    body_acceleration = 0 -9.81
  []
[]
[BoundaryConditions]
  [temperature]
    type = dirichlet
    boundary = all
    field = temperature
    value = 500
  []
  [fix_r]
    type = dirichlet
    boundary = left
    field = radial_displacement
    value = 0
  []
  [fix_z]
    type = dirichlet
    boundary = bottom
    field = axial_displacement
    value = 0
  []
[]
[Executioner]
  type = steady
  load_steps = 4
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-12
  relative_tolerance = 1e-12
  step_tolerance = 1e-12
  maximum_iterations = 40
[]
[Outputs]
  console = false
  csv = steady_gravity_quad8_summary.csv
  exodus = steady_gravity_quad8_results.e
[]
