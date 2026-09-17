[Case]
  version = 3
  problem = steady
  geometry = generalized_plane_strain
[]
[Mesh]
  type = exodus
  file = distorted.e
[]
[Materials]
  [solid]
    [thermal]
      function = linear_temperature_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
      reference_temperature = 300
      conductivity_temperature_coefficient = 0.01
      specific_heat_temperature_coefficient = 0.1
    []
    [elasticity]
      function = linear_temperature_isotropic
      young_modulus = 1e6
      poisson_ratio = 0.25
      reference_temperature = 300
      young_modulus_temperature_coefficient = 1000
      poisson_ratio_temperature_coefficient = 0
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
    u3 = 0.001
    rotation_x = 0.002
    rotation_y = -0.003
  []
[]
[BoundaryConditions]
  [bottom_x]
    type = dirichlet
    boundary = bottom
    field = displacement_x
    value = -0.00005
  []
  [top_x]
    type = dirichlet
    boundary = top
    field = displacement_x
    value = 0.00005
  []
  [middle_x]
    type = dirichlet
    boundary = middle
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
  [hot]
    type = dirichlet
    boundary = top
    field = temperature
    value = 400
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
  exodus = shear_distorted_results.e
[]
