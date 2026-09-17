[Case]
  version = 3
  problem = steady
  geometry = generalized_plane_strain
[]
[Mesh]
  type = exodus
  file = convergence_2.e
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
    volumetric_heat_source = 1000
  []
[]
[GeneralizedPlaneStrain]
  [solid]
    blocks = solid
    initial_thickness = 0.1
    u3 = 0
    rotation_x = 0
    rotation_y = 0
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
  [cold]
    type = dirichlet
    boundary = bottom
    field = temperature
    value = 300
  []
  [top_convection]
    type = convection
    boundary = top
    heat_transfer_coefficient = 10
    ambient_temperature = 350
  []
[]
[Executioner]
  type = steady
  load_steps = 1
[]
[Solver]
  absolute_tolerance = 1e-9
  relative_tolerance = 1e-11
  step_tolerance = 1e-14
  maximum_iterations = 30
  linear_solver = direct
  direct_factorization = mumps
[]
[Outputs]
  console = true
  exodus = mesh_convergence_2_results.e
[]
