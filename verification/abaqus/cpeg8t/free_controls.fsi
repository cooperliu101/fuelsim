[Case]
  version = 3
  problem = steady
  geometry = generalized_plane_strain
[]
[Mesh]
  type = exodus
  file = rectangle.e
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
      young_modulus = 1e6
      poisson_ratio = 0.25
    []
    [eigenstrains]
      [expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 300
      []
    []
  []
[]
[Regions]
  [solid]
    element = cpeg8t
    block = solid
    material = solid
    strain = small
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]
[GeneralizedPlaneStrain]
  [solid]
    blocks = solid
    initial_thickness = 0.1
  []
[]
[BoundaryConditions]
  [fix_x]
    type = dirichlet
    boundary = field
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = field
    field = displacement_y
    value = 0
  []
  [temperature1]
    type = dirichlet
    boundary = corner1
    field = temperature
    value = 390
  []
  [temperature2]
    type = dirichlet
    boundary = corner2
    field = temperature
    value = 410
  []
  [temperature3]
    type = dirichlet
    boundary = corner3
    field = temperature
    value = 420
  []
  [temperature4]
    type = dirichlet
    boundary = corner4
    field = temperature
    value = 400
  []
[]
[Executioner]
  type = steady
  load_steps = 1
[]
[Solver]
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-11
  step_tolerance = 1e-14
  maximum_iterations = 20
  linear_solver = direct
  direct_factorization = mumps
[]
[Outputs]
  console = true
  exodus = free_controls_results.e
[]
