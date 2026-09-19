[Case]
  version = 3
  physics = thermal
  problem = steady
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = pellet.e
[]
[Materials]
  [fuel]
    [thermal]
      function = constant_thermophysical
      conductivity = 3
      density = 10000
      specific_heat = 300
    []
  []
[]
[Regions]
  [pellet]
    block = pellet
    element = dc3d8
    material = fuel
    initial_temperature = 600
    volumetric_heat_source = 200000000
    pellet_response = exact
  []
[]
[BoundaryConditions]
  [surface_0]
    type = dirichlet
    boundary = surface_0
    field = temperature
    value = 664.6446609406726
  []
  [surface_1]
    type = dirichlet
    boundary = surface_1
    field = temperature
    value = 680
  []
  [surface_2]
    type = dirichlet
    boundary = surface_2
    field = temperature
    value = 707.0710678118655
  []
  [surface_3]
    type = dirichlet
    boundary = surface_3
    field = temperature
    value = 670
  []
  [surface_4]
    type = dirichlet
    boundary = surface_4
    field = temperature
    value = 700
  []
  [surface_5]
    type = dirichlet
    boundary = surface_5
    field = temperature
    value = 730
  []
  [surface_6]
    type = dirichlet
    boundary = surface_6
    field = temperature
    value = 692.9289321881345
  []
  [surface_7]
    type = dirichlet
    boundary = surface_7
    field = temperature
    value = 720
  []
  [surface_8]
    type = dirichlet
    boundary = surface_8
    field = temperature
    value = 735.3553390593274
  []
  [surface_9]
    type = dirichlet
    boundary = surface_9
    field = temperature
    value = 714.6446609406726
  []
  [surface_10]
    type = dirichlet
    boundary = surface_10
    field = temperature
    value = 730
  []
  [surface_11]
    type = dirichlet
    boundary = surface_11
    field = temperature
    value = 757.0710678118655
  []
  [surface_12]
    type = dirichlet
    boundary = surface_12
    field = temperature
    value = 720
  []
  [surface_14]
    type = dirichlet
    boundary = surface_14
    field = temperature
    value = 780
  []
  [surface_15]
    type = dirichlet
    boundary = surface_15
    field = temperature
    value = 742.9289321881345
  []
  [surface_16]
    type = dirichlet
    boundary = surface_16
    field = temperature
    value = 770
  []
  [surface_17]
    type = dirichlet
    boundary = surface_17
    field = temperature
    value = 785.3553390593274
  []
  [surface_18]
    type = dirichlet
    boundary = surface_18
    field = temperature
    value = 764.6446609406726
  []
  [surface_19]
    type = dirichlet
    boundary = surface_19
    field = temperature
    value = 780
  []
  [surface_20]
    type = dirichlet
    boundary = surface_20
    field = temperature
    value = 807.0710678118655
  []
  [surface_21]
    type = dirichlet
    boundary = surface_21
    field = temperature
    value = 770
  []
  [surface_22]
    type = dirichlet
    boundary = surface_22
    field = temperature
    value = 800
  []
  [surface_23]
    type = dirichlet
    boundary = surface_23
    field = temperature
    value = 830
  []
  [surface_24]
    type = dirichlet
    boundary = surface_24
    field = temperature
    value = 792.9289321881345
  []
  [surface_25]
    type = dirichlet
    boundary = surface_25
    field = temperature
    value = 820
  []
  [surface_26]
    type = dirichlet
    boundary = surface_26
    field = temperature
    value = 835.3553390593274
  []
[]
[Executioner]
  type = steady
  load_steps = 1
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 30
[]
[Outputs]
  console = true
  exodus = pellet_mixed_exact_results.e
[]
